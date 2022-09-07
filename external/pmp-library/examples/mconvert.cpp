// Copyright 2011-2021 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include <pmp/SurfaceMesh.h>
#include <unistd.h>

using namespace pmp;

void usage_and_exit()
{
    std::cerr << "Usage:\nmconvert [-b] -i <input> -o <output>\n\nOptions\n"
              << " -b:  write binary format\n"
              << "\n";
    exit(1);
}

int main(int argc, char** argv)
{
    bool binary = false;
    const char* input = nullptr;
    const char* output = nullptr;

    // parse command line parameters
    int c;
    while ((c = getopt(argc, argv, "bi:o:")) != -1)
    {
        switch (c)
        {
            case 'b':
                binary = true;
                break;

            case 'i':
                input = optarg;
                break;

            case 'o':
                output = optarg;
                break;

            default:
                usage_and_exit();
        }
    }

    // we need input and output mesh
    if (!input || !output)
    {
        usage_and_exit();
    }

    // load input mesh
    SurfaceMesh mesh;
    try
    {
        mesh.read(input);
    }
    catch (const IOException& e)
    {
        std::cerr << "Failed to read mesh: " << e.what() << std::endl;
        exit(1);
    }

    // write output mesh
    IOFlags flags;
    flags.use_binary = binary;
    try
    {
        mesh.write(output, flags);
    }
    catch (const IOException& e)
    {
        std::cerr << "Failed to write mesh: " << e.what() << std::endl;
        exit(1);
    }

    exit(0);
}
