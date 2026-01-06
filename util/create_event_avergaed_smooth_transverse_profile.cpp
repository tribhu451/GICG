#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

struct Cell {
    double mult = 0.0;
    double mult_a = 0.0;
    double mult_b = 0.0;
    double mult_coll = 0.0;
    double npart_a = 0.0;
    double npart_b = 0.0;
    double ncoll = 0.0;
};

int main() {

    const int nEvents = 10;  // number of events to average

    int nx = -1, ny = -1;
    double dx = 0.0, dy = 0.0;

    std::vector<std::vector<Cell>> avg;

    for (int ievt = 0; ievt < nEvents; ++ievt) {

        std::stringstream filename;
        filename << "../output/mc_glauber_single_event_transverse_profile_for_rapidity_extension_"
                 << ievt << ".dat";

        std::ifstream infile(filename.str());
        if (!infile.is_open()) {
            std::cerr << "Cannot open " << filename.str() << std::endl;
            continue;
        }

        std::string line;

        // ---------- Read header ----------
        std::getline(infile, line);
        if (line[0] != '#') {
            std::cerr << "Header missing in " << filename.str() << std::endl;
            return 1;
        }

        if (ievt == 0) {
            std::istringstream h(line);
            std::string tmp;

            h >> tmp >> tmp >> tmp
              >> tmp >> tmp
              >> tmp >> nx
              >> tmp >> ny
              >> tmp >> tmp
              >> tmp >> dx
              >> tmp >> dy;

            avg.resize(nx, std::vector<Cell>(ny));
        }

        // ---------- Read transverse cells ----------
        while (std::getline(infile, line)) {

            if (line.empty() || line[0] == '#') continue;

            std::istringstream iss(line);

            int ieta;
            double x, y;
            Cell c;
            double dummy;

            iss >> ieta >> x >> y
                >> c.mult >> c.mult_a >> c.mult_b >> c.mult_coll
                >> c.npart_a >> c.npart_b >> c.ncoll
                >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;

            int ix = static_cast<int>((x + 0.5 * nx * dx) / dx);
            int iy = static_cast<int>((y + 0.5 * ny * dy) / dy);

            if (ix < 0 || ix >= nx || iy < 0 || iy >= ny) continue;

            avg[ix][iy].mult      += c.mult;
            avg[ix][iy].mult_a    += c.mult_a;
            avg[ix][iy].mult_b    += c.mult_b;
            avg[ix][iy].mult_coll += c.mult_coll;
            avg[ix][iy].npart_a   += c.npart_a;
            avg[ix][iy].npart_b   += c.npart_b;
            avg[ix][iy].ncoll     += c.ncoll;
        }

        infile.close();
    }

    // ---------- Event average ----------
    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            avg[ix][iy].mult      /= nEvents;
            avg[ix][iy].mult_a    /= nEvents;
            avg[ix][iy].mult_b    /= nEvents;
            avg[ix][iy].mult_coll /= nEvents;
            avg[ix][iy].npart_a   /= nEvents;
            avg[ix][iy].npart_b   /= nEvents;
            avg[ix][iy].ncoll     /= nEvents;
        }
    }

    // ---------- Write averaged transverse profile ----------
    std::ofstream out("event_averaged_transverse_profile_eta0.dat");

    out << "# neta 1 nx " << nx << " ny " << ny
        << " dx " << dx << " dy " << dy << "\n";

    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {

            double x = -0.5 * (nx-1) * dx + ix * dx;
            double y = -0.5 * (ny-1) * dy + iy * dy;

            out << "0 " << x << " " << y << " "
                << avg[ix][iy].mult << " "
                << avg[ix][iy].mult_a << " "
                << avg[ix][iy].mult_b << " "
                << avg[ix][iy].mult_coll << " "
                << avg[ix][iy].npart_a << " "
                << avg[ix][iy].npart_b << " "
                << avg[ix][iy].ncoll << "\n";
        }
    }

    out.close();

    std::cout << "Event-averaged transverse profile at eta=0 written.\n";
    return 0;
}

