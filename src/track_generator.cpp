#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <cstdlib>

constexpr int N_DETECTORS = 8;
constexpr int N_EVENTS = 10000;
constexpr float BEAM_SIGMA = 1.0f;          // cm, beam spot size in X and Y
constexpr float BEAM_DIVERGENCE = 0.001f;   // rad, angular spread (~1 mrad)
constexpr float BEAM_CENTER_X = 0.0f;       // cm
constexpr float BEAM_CENTER_Y = 0.0f;       // cm
constexpr float DETECTOR_WIDTH = 20.0f;     // cm, full width (+-10 cm)
constexpr float LAYER_SPACING = 10.0f;      // cm, nominal z spacing
constexpr float MEASUREMENT_NOISE = 0.003f; // cm, 30 um spatial resolution
constexpr int MIN_HITS = 4;

// ---------------------------------------------------------------------
// Minimal 3D vector / rotation matrix helpers
// ---------------------------------------------------------------------
struct Vec3
{
    float x = 0, y = 0, z = 0;
    Vec3 operator+(const Vec3 &o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3 &o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(float s) const { return {x * s, y * s, z * s}; }
    float dot(const Vec3 &o) const { return x * o.x + y * o.y + z * o.z; }
};

struct Mat3
{
    float m[3][3] = {};
    Vec3 apply(const Vec3 &v) const
    {
        return {
            m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z};
    }
};

static Mat3 matmul(const Mat3 &A, const Mat3 &B)
{
    Mat3 C;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            C.m[i][j] = 0.0f;
            for (int k = 0; k < 3; k++)
                C.m[i][j] += A.m[i][k] * B.m[k][j];
        }
    return C;
}

// Full rotation R = Rz * Ry * Rx (Rx applied first, then Ry, then Rz)
static Mat3 rotation_matrix(float rx, float ry, float rz)
{
    float cx = std::cos(rx), sx = std::sin(rx);
    float cy = std::cos(ry), sy = std::sin(ry);
    float cz = std::cos(rz), sz = std::sin(rz);

    Mat3 Rx;
    Rx.m[0][0] = 1;
    Rx.m[0][1] = 0;
    Rx.m[0][2] = 0;
    Rx.m[1][0] = 0;
    Rx.m[1][1] = cx;
    Rx.m[1][2] = -sx;
    Rx.m[2][0] = 0;
    Rx.m[2][1] = sx;
    Rx.m[2][2] = cx;

    Mat3 Ry;
    Ry.m[0][0] = cy;
    Ry.m[0][1] = 0;
    Ry.m[0][2] = sy;
    Ry.m[1][0] = 0;
    Ry.m[1][1] = 1;
    Ry.m[1][2] = 0;
    Ry.m[2][0] = -sy;
    Ry.m[2][1] = 0;
    Ry.m[2][2] = cy;

    Mat3 Rz;
    Rz.m[0][0] = cz;
    Rz.m[0][1] = -sz;
    Rz.m[0][2] = 0;
    Rz.m[1][0] = sz;
    Rz.m[1][1] = cz;
    Rz.m[1][2] = 0;
    Rz.m[2][0] = 0;
    Rz.m[2][1] = 0;
    Rz.m[2][2] = 1;

    return matmul(matmul(Rz, Ry), Rx);
}

// ---------------------------------------------------------------------
// Detector: nominal geometry + 6 DOF alignment
// ---------------------------------------------------------------------
enum class ReadoutAxis
{
    X,
    Y
};

struct Detector
{
    // --- nominal (design) geometry ---
    Vec3 nominal_position;    // usually {0, 0, z_layer}
    ReadoutAxis readout_axis; // which local coordinate this plane measures
    float width = DETECTOR_WIDTH;

    // --- 6 DOF alignment on top of nominal geometry ---
    float dx = 0, dy = 0, dz = 0; // translations [cm]
    float rx = 0, ry = 0, rz = 0; // rotations around X,Y,Z [rad]

    // --- derived quantities, computed once by finalize() ---
    Vec3 center; // actual plane center
    Vec3 normal; // actual plane normal
    Vec3 u_axis; // actual in-plane readout direction

    void finalize()
    {
        center = nominal_position + Vec3{dx, dy, dz};
        Mat3 R = rotation_matrix(rx, ry, rz);
        normal = R.apply(Vec3{0, 0, 1});
        Vec3 nominal_u = (readout_axis == ReadoutAxis::X) ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
        u_axis = R.apply(nominal_u);
    }
};

// ---------------------------------------------------------------------
// Track: straight line, origin at z=0, direction ~ (slope_x, slope_y, 1)
// ---------------------------------------------------------------------
struct Track
{
    Vec3 origin;
    Vec3 direction;
};

// Intersect the track with a detector plane and return the local readout
// coordinate (signed distance from plane center along u_axis).
// Returns false only if the track runs parallel to the plane.
static bool project_track(const Track &trk, const Detector &det, float &local_u)
{
    float denom = trk.direction.dot(det.normal);
    if (std::fabs(denom) < 1e-9f)
        return false;

    float t = (det.center - trk.origin).dot(det.normal) / denom;
    Vec3 hit = trk.origin + trk.direction * t;
    local_u = (hit - det.center).dot(det.u_axis);
    return true;
}

// ---------------------------------------------------------------------
// Alignment file: one entry per misaligned layer
//   layer,dx[cm],dy[cm],dz[cm],rx[mrad],ry[mrad],rz[mrad]
// Lines starting with # (after stripping) are comments; blank lines and
// a non-numeric header row are skipped automatically. Layers not listed
// keep the ideal (zero) alignment.
// ---------------------------------------------------------------------
struct AlignmentEntry
{
    int layer;
    float dx, dy, dz; // cm
    float rx, ry, rz; // rad (converted from mrad on read)
};

static std::vector<AlignmentEntry> read_alignment_file(const std::string &path)
{
    std::vector<AlignmentEntry> entries;
    if (path.empty())
        return entries;

    std::ifstream f(path);
    if (!f.is_open())
    {
        std::cerr << "WARNING: alignment file '" << path
                  << "' not found, using ideal geometry" << std::endl;
        return entries;
    }

    std::string line;
    int line_no = 0;
    while (std::getline(f, line))
    {
        line_no++;
        auto hash = line.find('#');
        if (hash != std::string::npos)
            line = line.substr(0, hash);
        // trim leading/trailing whitespace
        size_t start = line.find_first_not_of(" \t\r\n");
        if (start == std::string::npos)
            continue; // blank/comment-only line
        size_t end = line.find_last_not_of(" \t\r\n");
        line = line.substr(start, end - start + 1);

        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> fields;
        while (std::getline(ss, token, ','))
            fields.push_back(token);
        if (fields.size() != 7)
        {
            std::cerr << "WARNING: line " << line_no << " in " << path
                      << " has " << fields.size() << " fields (expected 7), skipping"
                      << std::endl;
            continue;
        }

        try
        {
            AlignmentEntry e;
            e.layer = std::stoi(fields[0]);
            e.dx = std::stof(fields[1]);
            e.dy = std::stof(fields[2]);
            e.dz = std::stof(fields[3]);
            e.rx = std::stof(fields[4]) * 0.001f; // mrad -> rad
            e.ry = std::stof(fields[5]) * 0.001f;
            e.rz = std::stof(fields[6]) * 0.001f;
            entries.push_back(e);
        }
        catch (const std::exception &)
        {
            // Non-numeric row (e.g. header like "layer,dx,dy,..."): skip silently.
            continue;
        }
    }
    return entries;
}

// ---------------------------------------------------------------------
// Detector array definition
// ---------------------------------------------------------------------
static std::array<Detector, N_DETECTORS> build_detectors(const std::vector<AlignmentEntry> &alignment)
{
    std::array<Detector, N_DETECTORS> dets;
    for (int i = 0; i < N_DETECTORS; i++)
    {
        Detector d;
        d.nominal_position = Vec3{0.0f, 0.0f, i * LAYER_SPACING};
        d.readout_axis = (i % 2 == 0) ? ReadoutAxis::X : ReadoutAxis::Y;
        d.width = DETECTOR_WIDTH;
        dets[i] = d; // 6 DOF default to zero (ideal); applied below, finalized after
    }

    for (const auto &e : alignment)
    {
        if (e.layer < 0 || e.layer >= N_DETECTORS)
        {
            std::cerr << "WARNING: alignment entry for layer " << e.layer
                      << " out of range [0," << N_DETECTORS - 1 << "], ignoring" << std::endl;
            continue;
        }
        dets[e.layer].dx = e.dx;
        dets[e.layer].dy = e.dy;
        dets[e.layer].dz = e.dz;
        dets[e.layer].rx = e.rx;
        dets[e.layer].ry = e.ry;
        dets[e.layer].rz = e.rz;
    }
    for (auto &d : dets)
        d.finalize();
    return dets;
}

// ---------------------------------------------------------------------
int main(int argc, char *argv[])
{

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <n_events> [<beam_sigma> [<out_file> [<divergence> [<align_file>]]]]" << std::endl;
        return 1;
    }

    int n_events = (argc > 1) ? std::atoi(argv[1]) : N_EVENTS;
    float beam_sigma = (argc > 2) ? std::atof(argv[2]) : BEAM_SIGMA;
    std::string out_file = (argc > 3) ? argv[3] : "hits_generated.csv";
    float divergence = (argc > 4) ? std::atof(argv[4]) : BEAM_DIVERGENCE;
    std::string align_file = (argc > 5) ? argv[5] : "";

    std::cout << "Generating " << n_events << " events "
              << "(beam sigma=" << beam_sigma << " cm, "
              << "divergence=" << divergence << " rad)" << std::endl;

    auto alignment = read_alignment_file(align_file);
    if (!align_file.empty() && !alignment.empty())
    {
        std::cout << "Applying " << alignment.size() << " alignment entr"
                  << (alignment.size() == 1 ? "y" : "ies")
                  << " from " << align_file << std::endl;
    }

    if (!align_file.empty() && !alignment.empty())
    {
        // Print all alignment entries
        std::cout << "\nAlignment Entries:" << std::endl;
        for (const auto &e : alignment)
        {
            std::cout << std::left << std::setw(20) << e.layer
                      << std::right << std::setw(15) << std::fixed
                      << std::setprecision(6) << e.dx
                      << std::right << std::setw(15) << std::fixed
                      << std::setprecision(6) << e.dy
                      << std::right << std::setw(15) << std::fixed
                      << std::setprecision(6) << e.dz
                      << std::right << std::setw(15) << std::fixed
                      << std::setprecision(6) << e.rx
                      << std::right << std::setw(15) << std::fixed
                      << std::setprecision(6) << e.ry
                      << std::right << std::setw(15) << std::fixed
                      << std::setprecision(6) << e.rz << std::endl;
        }
    }
    auto detectors = build_detectors(alignment);

    std::mt19937 gen(42);
    std::normal_distribution<float> beam_pos_dist(0.0f, beam_sigma);
    std::normal_distribution<float> beam_angle_dist(0.0f, divergence);
    std::normal_distribution<float> noise_dist(0.0f, MEASUREMENT_NOISE);

    std::ofstream csv(out_file);
    if (!csv.is_open())
    {
        std::cerr << "ERROR: Cannot open " << out_file << std::endl;
        return 1;
    }

    int valid_events = 0;
    for (int event = 0; event < n_events; event++)
    {
        Track trk;
        trk.origin = Vec3{
            BEAM_CENTER_X + beam_pos_dist(gen),
            BEAM_CENTER_Y + beam_pos_dist(gen),
            0.0f};
        trk.direction = Vec3{
            beam_angle_dist(gen),
            beam_angle_dist(gen),
            1.0f};

        std::array<float, N_DETECTORS> hits;
        for (int det = 0; det < N_DETECTORS; det++)
        {
            float local_u;
            bool ok = project_track(trk, detectors[det], local_u);
            if (!ok)
            {
                hits[det] = std::numeric_limits<float>::quiet_NaN();
                continue;
            }
            local_u += noise_dist(gen);

            if (std::fabs(local_u) > detectors[det].width / 2.0f)
            {
                hits[det] = std::numeric_limits<float>::quiet_NaN();
            }
            else
            {
                hits[det] = local_u;
            }
        }

        int n_hits = 0;
        for (int det = 0; det < N_DETECTORS; det++)
            if (!std::isnan(hits[det]))
                n_hits++;

        if (n_hits < MIN_HITS)
            continue;

        valid_events++;
        for (int det = 0; det < N_DETECTORS; det++)
        {
            if (det > 0)
                csv << ",";
            if (std::isnan(hits[det]))
                csv << "NaN";
            else
                csv << std::fixed << std::setprecision(6) << hits[det];
        }
        csv << "\n";
    }

    csv.close();
    std::cout << "Wrote " << valid_events << " valid events to " << out_file << std::endl;
    std::cout << "Acceptance: " << (100.0f * valid_events / n_events) << "%" << std::endl;
    return 0;
}
