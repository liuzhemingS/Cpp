#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace constants {
constexpr double kCoulomb = 8.9875517923e9;      // N m^2 / C^2
constexpr double eCharge = -1.602176634e-19;     // C (electron charge)
constexpr double eMass = 9.1093837015e-31;       // kg
constexpr double c = 299792458.0;                // m / s
constexpr double kBoltzmann = 1.380649e-23;      // J / K
constexpr double epsilon0 = 8.8541878128e-12;    // F / m
} // namespace constants

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& rhs) const { return {x + rhs.x, y + rhs.y, z + rhs.z}; }
    Vec3 operator-(const Vec3& rhs) const { return {x - rhs.x, y - rhs.y, z - rhs.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }

    Vec3& operator+=(const Vec3& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    Vec3& operator-=(const Vec3& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }
};

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    };
}

double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

double norm(const Vec3& v) {
    return std::sqrt(dot(v, v));
}

struct Electron {
    Vec3 position;          // m
    Vec3 properVelocity;    // gamma * v, m/s
};

struct SimulationConfig {
    std::size_t particleCount = 200;
    std::size_t steps = 25000;
    std::size_t sampleEvery = 250;
    double dt = 5.0e-14;                      // s
    double cloudSigma = 5.0e-6;              // m
    double temperature = 1200.0;             // K
    double softeningLength = 2.0e-10;        // m

    // Idealized Penning-like confinement field settings.
    double trapElectricGradient = 8.0e8;     // V / m^2
    double uniformBz = 4.0;                  // Tesla
};

Vec3 externalElectricField(const Vec3& r, const SimulationConfig& cfg) {
    // Quadrupole electrostatic trap potential phi = (k/2)(z^2 - (x^2+y^2)/2)
    // E = -grad(phi)
    return {
        cfg.trapElectricGradient * 0.5 * r.x,
        cfg.trapElectricGradient * 0.5 * r.y,
        -cfg.trapElectricGradient * r.z,
    };
}

Vec3 externalMagneticField(const Vec3&, const SimulationConfig& cfg) {
    return {0.0, 0.0, cfg.uniformBz};
}

double gammaFromProperVelocity(const Vec3& u) {
    const double u2 = dot(u, u);
    return std::sqrt(1.0 + u2 / (constants::c * constants::c));
}

std::vector<Electron> initializeElectrons(const SimulationConfig& cfg, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> posDist(0.0, cfg.cloudSigma);

    // Maxwell-Boltzmann scale for each component of velocity.
    const double sigmaV = std::sqrt(constants::kBoltzmann * cfg.temperature / constants::eMass);
    std::normal_distribution<double> velDist(0.0, sigmaV);

    std::vector<Electron> electrons;
    electrons.reserve(cfg.particleCount);

    for (std::size_t i = 0; i < cfg.particleCount; ++i) {
        Vec3 v{velDist(rng), velDist(rng), velDist(rng)};
        const double speed = norm(v);
        if (speed > 0.98 * constants::c) {
            v = v * ((0.98 * constants::c) / speed);
        }

        const double gamma = 1.0 / std::sqrt(1.0 - dot(v, v) / (constants::c * constants::c));
        electrons.push_back({
            {posDist(rng), posDist(rng), posDist(rng)},
            v * gamma,
        });
    }

    return electrons;
}

std::vector<Vec3> computePairwiseElectricField(const std::vector<Electron>& electrons,
                                               const SimulationConfig& cfg) {
    std::vector<Vec3> fields(electrons.size());

    // O(N^2) Coulomb field from all other electrons.
    for (std::size_t i = 0; i < electrons.size(); ++i) {
        for (std::size_t j = i + 1; j < electrons.size(); ++j) {
            const Vec3 dr = electrons[i].position - electrons[j].position;
            const double r2 = dot(dr, dr) + cfg.softeningLength * cfg.softeningLength;
            const double invR = 1.0 / std::sqrt(r2);
            const double invR3 = invR * invR * invR;

            // E_i due to j and E_j due to i from point charge model with softening.
            const Vec3 eOnI = dr * (constants::kCoulomb * constants::eCharge * invR3);
            fields[i] += eOnI;
            fields[j] -= eOnI;
        }
    }

    return fields;
}

double kineticEnergyJ(const Electron& e) {
    const double gamma = gammaFromProperVelocity(e.properVelocity);
    return (gamma - 1.0) * constants::eMass * constants::c * constants::c;
}

double pairPotentialEnergyJ(const std::vector<Electron>& electrons, const SimulationConfig& cfg) {
    double pe = 0.0;
    for (std::size_t i = 0; i < electrons.size(); ++i) {
        for (std::size_t j = i + 1; j < electrons.size(); ++j) {
            const Vec3 dr = electrons[i].position - electrons[j].position;
            const double r = std::sqrt(dot(dr, dr) + cfg.softeningLength * cfg.softeningLength);
            pe += constants::kCoulomb * constants::eCharge * constants::eCharge / r;
        }
    }
    return pe;
}

int main() {
    SimulationConfig cfg;
    const std::uint64_t seed = 123456789ULL;

    std::vector<Electron> electrons = initializeElectrons(cfg, seed);

    std::ofstream csv("electron_diagnostics.csv");
    csv << "step,time_s,mean_kinetic_eV,total_kinetic_eV,total_pair_potential_eV,mean_radius_um\n";

    for (std::size_t step = 0; step <= cfg.steps; ++step) {
        const auto pairFields = computePairwiseElectricField(electrons, cfg);

        for (std::size_t i = 0; i < electrons.size(); ++i) {
            const Vec3 E = pairFields[i] + externalElectricField(electrons[i].position, cfg);
            const Vec3 B = externalMagneticField(electrons[i].position, cfg);

            // Relativistic Boris push in terms of proper velocity u = gamma*v.
            const double qOverM = constants::eCharge / constants::eMass;
            const Vec3 uMinus = electrons[i].properVelocity + E * (qOverM * cfg.dt * 0.5);
            const double gammaMinus = gammaFromProperVelocity(uMinus);

            const Vec3 t = B * (qOverM * cfg.dt / (2.0 * gammaMinus));
            const double t2 = dot(t, t);
            const Vec3 s = t * (2.0 / (1.0 + t2));

            const Vec3 uPrime = uMinus + cross(uMinus, t);
            const Vec3 uPlus = uMinus + cross(uPrime, s);
            const Vec3 uNew = uPlus + E * (qOverM * cfg.dt * 0.5);

            const double gammaNew = gammaFromProperVelocity(uNew);
            const Vec3 vNew = uNew / gammaNew;

            electrons[i].properVelocity = uNew;
            electrons[i].position += vNew * cfg.dt;
        }

        if (step % cfg.sampleEvery == 0) {
            double totalKE = 0.0;
            double meanRadius = 0.0;
            for (const auto& e : electrons) {
                totalKE += kineticEnergyJ(e);
                meanRadius += norm(e.position);
            }
            meanRadius /= static_cast<double>(electrons.size());

            const double totalPE = pairPotentialEnergyJ(electrons, cfg);
            const double ev = std::abs(constants::eCharge);

            csv << step << ','
                << std::setprecision(12) << step * cfg.dt << ','
                << (totalKE / electrons.size()) / ev << ','
                << totalKE / ev << ','
                << totalPE / ev << ','
                << meanRadius * 1e6 << '\n';

            std::cout << "step=" << step
                      << " meanKE[eV]=" << (totalKE / electrons.size()) / ev
                      << " meanRadius[um]=" << meanRadius * 1e6
                      << '\n';
        }
    }

    std::cout << "Simulation finished. Diagnostics saved to electron_diagnostics.csv\n";
    std::cout << "Model notes: classical electrodynamics with relativistic motion, Coulomb interactions,\n"
              << "external trap fields, and softening length to avoid singular close-encounters.\n";

    return 0;
}
