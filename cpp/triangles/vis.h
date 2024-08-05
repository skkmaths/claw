int fileid = 0;
double length_face(Face& face)
{
    Node* p0 = face.nodes[0];
    Node* p1 = face.nodes[1];
    return std::sqrt( std::pow (p0->x-p1->x,2) + std::pow (p0->y - p1->y,2) );
}
double minfacelength(const Mesh& mesh)
{
    double h = 1e20;
    for (auto &face : mesh.faces)
    {
        h = std::min( h, face.length);
    }
    return h;
}
// For file name
void createDirectory(const std::string& dirname) 
{
    struct stat info;
    if(stat(dirname.c_str(), &info) != 0 || !(info.st_mode & S_IFDIR)) {
        if(mkdir(dirname.c_str(), 0777) != 0) {
            std::cerr << "Error creating directory " << dirname << std::endl;
            exit(EXIT_FAILURE);
        }
        std::cout << "Directory " << dirname << " is created" << std::endl;
    }
}
std::string getFilename(const std::string& basename, int id) {
    std::ostringstream oss;
    oss << basename << "_" << std::setfill('0') << std::setw(4) << id << ".vtk";
    return oss.str();
}
// Function to write initial condition to a VTK file
void savesol(const Mesh& mesh, std::vector<double>& solution, const double& t) {
    std::string dirname = "sol";
    createDirectory(dirname);
    if(fileid == 0) {
        std::cout << "The directory \"sol\" is going to be formatted!" << std::endl;
        std::string pattern = "./sol/*";
        // Remove existing files
        system(("rm -f " + pattern).c_str());
    }
    std::string filename = getFilename("sol/sol", fileid);
    std::ofstream file(filename);
    file << "# vtk DataFile Version 3.0\n";
    file << "Advection Solution at time " << t << "\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n";

    // Write points
    // Warning! make sure that id s are index from 0 to mesh.nodes.size()
    file << "POINTS " << mesh.nodes.size() << " float\n";
    for ( unsigned int id = 0; id < mesh.nodes.size(); id++)
    {
        file <<mesh.nodeMap[id]->x<<" " <<mesh.nodeMap[id]->y <<" "<< mesh.nodeMap[id]->z << "\n";
    }

    // Write cells
    file << "CELLS " << mesh.cells.size() << " " << 4 * mesh.cells.size() << "\n";
    for (const auto &cell : mesh.cells) {
        file << "3 " << cell.nodes[0]->id << " " << cell.nodes[1]->id << " " << cell.nodes[2]->id << "\n";
    }

    // Write cell types
    file << "CELL_TYPES " << mesh.cells.size() << "\n";
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        file << "5\n"; // VTK_TRIANGLE
    }

    // Write field data for time
    file << "FIELD FieldData 1\n";
    file << "TIME 1 1 float\n";
    file << t << "\n";

    // Write cell data
    file << "CELL_DATA " << mesh.cells.size() << "\n";
    file << "SCALARS sol float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto &cell : mesh.cells) {
        file << std::fixed << std::setprecision(12) << solution[cell.id] << "\n";
    }

    file.close();
    fileid++;
}