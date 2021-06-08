#include "MarchingCubes.h"
#include "MarchingCubes.hh"

namespace Ultraliser
{

MarchingCubes::MarchingCubes(Volume* volume,
                             const uint8_t _isoValue)
    : _volume(volume)
    , _isoValue(_isoValue)
{
    /// EMPTY CONSTRUCTOR
}






double mc_isoValue_interpolation(double _isoValue, double f1, double f2,
                                double x1, double x2)
{
    if(f2==f1)
        return (x2+x1)/2;

    return (x2-x1)*(_isoValue-f1)/(f2-f1) + x1;
}


size_t mc_add_vertex3(double x1, double y1, double z1, double c2,
    int axis, double f1, double f2, double _isoValue, Vertices& vertices)
{
   // size_t vertex_index = vertices->size() / 3;

    if(axis == 0)
    {
        double x = mc_isoValue_interpolation(_isoValue, f1, f2, x1, c2);
        vertices.push_back(Vector3f(x, y1, z1));
//        vertices->push_back(y1);
//        vertices->push_back(z1);
 //       return vertex_index;
    }
    if(axis == 1)
    {
        double y = mc_isoValue_interpolation(_isoValue, f1, f2, y1, c2);
        vertices.push_back(Vector3f(x1, y, z1));
//        vertices->push_back(y);
//        vertices->push_back(z1);
   //     return vertex_index;
    }
    if(axis == 2)
    {
        double z = mc_isoValue_interpolation(_isoValue, f1, f2, z1, c2);
        vertices.push_back(Vector3f(x1, y1, z));
//        vertices->push_back(y1);
//        vertices->push_back(z);
  //      return vertex_index;
    }

    // This should not happen.
    return -1;
}



size_t mc_add_vertex(double x1, double y1, double z1, double c2,
    int axis, double f1, double f2, double _isoValue, std::vector<double>* vertices)
{
    size_t vertex_index = vertices->size() / 3;

    if(axis == 0)
    {
        double x = mc_isoValue_interpolation(_isoValue, f1, f2, x1, c2);
        vertices->push_back(x);
        vertices->push_back(y1);
        vertices->push_back(z1);
        return vertex_index;
    }
    if(axis == 1)
    {
        double y = mc_isoValue_interpolation(_isoValue, f1, f2, y1, c2);
        vertices->push_back(x1);
        vertices->push_back(y);
        vertices->push_back(z1);
        return vertex_index;
    }
    if(axis == 2)
    {
        double z = mc_isoValue_interpolation(_isoValue, f1, f2, z1, c2);
        vertices->push_back(x1);
        vertices->push_back(y1);
        vertices->push_back(z);
        return vertex_index;
    }

    // This should not happen.
    return -1;
}

//template<typename vector3, typename formula>
//void marching_cubes(const vector3& lower, const vector3& upper,
//    int numx, int numy, int numz, formula f, double _isoValue,
//    std::vector<double>& vertices, std::vector<typename vector3::size_type>& polygons)


Mesh* MarchingCubes::marching_cubes()

{
    /// f is the function and _isoValue
    /// pmin and pmax
    /// const vector3& lower, const vector3& upper
    ///
    /// volume dimensions
    /// int numx, int numy, int numz,
    ///
    //// std::vector<double>& vertices, std::vector<typename vector3::size_type>& polygons
    ///
    // Build the mesh
    Vertices verticesA;
    Triangles trianglesA;

    std::vector<size_t> polygons;
    std::vector<double> vertices;

   //  template<typename vector3, typename formula>
   // using coord_type = typename vector3::value_type;
   // using size_type = typename vector3::size_type;



//    // Some initial checks
//    if(numx < 2 || numy < 2 || numz < 2)
//        return;

//    if(!std::equal(std::begin(lower), std::end(lower), std::begin(upper),
//                   [](double a, double b)->bool {return a <= b;}))
//        return;

//    // numx, numy and numz are the numbers of evaluations in each direction
//    --numx; --numy; --numz;

//    coord_type dx = (upper[0] - lower[0]) / static_cast<coord_type>(numx);
//    coord_type dy = (upper[1] - lower[1]) / static_cast<coord_type>(numy);
//    coord_type dz = (upper[2] - lower[2]) / static_cast<coord_type>(numz);

    // Delta between the voxels
    uint64_t dx = 1;
    uint64_t dy = 1;
    uint64_t dz = 1;



    int64_t extraVoxels = 5;
    int64_t minValue = -1 * extraVoxels;
    int64_t maxValue = extraVoxels;

    uint64_t numx = _volume->getWidth() + maxValue;
    uint64_t numy = _volume->getHeight() + maxValue;
    uint64_t numz = _volume->getDepth() + maxValue;

    uint64_t size = ((numx + extraVoxels) * (numy + extraVoxels) * (numz + extraVoxels) * 3);

    size_t* shared_indices = new size_t[size];
    const int z3 = numz*3;
    const int yz3 = numy*z3;

    for (int i = minValue; i < _volume->getWidth() + 5; ++i)
    {
        int64_t x = dx * i;
        int64_t x_dx = dx * (i + 1);

        for (int j = minValue; j < _volume->getHeight() + 5; ++j)
        {
            int64_t y = dy*j;
            int64_t y_dy = dy * (j + 1);

            for (int k = minValue; k < _volume->getDepth() + 5; ++k)
            {
                int64_t z = dz*k;
                int64_t z_dz = dz*(k+1);

                double v[8];
                v[0] = _volume->getValue(i, j, k);
                v[1] = _volume->getValue(i + 1, j, k);
                v[2] = _volume->getValue(i + 1, j + 1, k);
                v[3] = _volume->getValue(i, j + 1, k);
                v[4] = _volume->getValue(i, j, k + 1);
                v[5] = _volume->getValue(i + 1, j, k + 1);
                v[6] = _volume->getValue(i + 1, j + 1, k + 1);
                v[7] = _volume->getValue(i, j + 1, k + 1);

                unsigned int cubeindex = 0;
                for(int m=0; m<8; ++m)
                    if(v[m] <= _isoValue)
                        cubeindex |= 1<<m;

                // Generate vertices AVOIDING DUPLICATES.

                int edges = edge_table[cubeindex];
                std::array<size_t, 12> indices;
                if(edges & 0x040)
                {
                    indices[6] = vertices.size() / 3;
                    int64_t idx = i*yz3 + j*z3 + k*3 + 0;
                    if (idx < 0) continue;
                    shared_indices[idx] = indices[6];
                    mc_add_vertex(x_dx, y_dy, z_dz, x, 0, v[6], v[7], _isoValue, &vertices);
                }
                if(edges & 0x020)
                {

                    indices[5] = vertices.size() / 3;

                    int64_t idx = i*yz3 + j*z3 + k*3 + 1;
                    if (idx < 0) continue;
                    shared_indices[idx] = indices[5];
                    mc_add_vertex(x_dx, y, z_dz, y_dy, 1, v[5], v[6], _isoValue, &vertices);
                }
                if(edges & 0x400)
                {


                    indices[10] = vertices.size() / 3;

                    int64_t idx = i*yz3 + j*z3 + k*3 + 2;
                    if (idx < 0) continue;
                    shared_indices[idx] = indices[10];
                    mc_add_vertex(x_dx, y+dx, z, z_dz, 2, v[2], v[6], _isoValue, &vertices);
                }

                if(edges & 0x001)
                {
                    if(j == 0 || k == 0)
                    {

                        indices[0] = vertices.size() / 3;
                        mc_add_vertex(x, y, z, x_dx, 0, v[0], v[1], _isoValue, &vertices);
                    }
                    else
                    {
                        int64_t idx = i*yz3 + (j-1)*z3 + (k-1)*3 + 0;
                        if (idx < 0) continue;
                        indices[0] = shared_indices[idx];
                    }
                }
                if(edges & 0x002)
                {
                    if(k == 0)
                    {
                        indices[1] = vertices.size() / 3;
                        mc_add_vertex(x_dx, y, z, y_dy, 1, v[1], v[2], _isoValue, &vertices);
                    }
                    else
                    {
                        int64_t idx = i*yz3 + j*z3 + (k-1)*3 + 1;
                        if (idx < 0) continue;
                        indices[1] = shared_indices[idx];
                    }
                }
                if(edges & 0x004)
                {
                    if(k == 0)
                    {
                        indices[2] = vertices.size() / 3;
                        mc_add_vertex(x_dx, y_dy, z, x, 0, v[2], v[3], _isoValue, &vertices);
                    }
                    else
                    {
                        int64_t idx = i*yz3 + j*z3 + (k-1)*3 + 0;
                        if (idx < 0) continue;
                        indices[2] = shared_indices[idx];
                    }
                }
                if(edges & 0x008)
                {
                    if(i == 0 || k == 0)
                    {
                        indices[3] = vertices.size() / 3;
                        mc_add_vertex(x, y_dy, z, y, 1, v[3], v[0], _isoValue, &vertices);
                    }
                    else
                    {
                        int64_t idx = (i-1)*yz3 + j*z3 + (k-1)*3 + 1;
                        if (idx < 0) continue;
                        indices[3] = shared_indices[idx];
                    }
                }
                if(edges & 0x010)
                {
                    if(j == 0)
                    {
                        indices[4] = vertices.size() / 3;
                        mc_add_vertex(x, y, z_dz, x_dx, 0, v[4], v[5], _isoValue, &vertices);
                    }
                    else
                    {
                        int64_t idx = i*yz3 + (j-1)*z3 + k*3 + 0;
                        if (idx < 0) continue;
                        indices[4] = shared_indices[idx];
                    }
                }
                if(edges & 0x080)
                {
                    if(i == 0)
                    {
                        indices[7] = vertices.size() / 3;
                        mc_add_vertex(x, y_dy, z_dz, y, 1, v[7], v[4], _isoValue, &vertices);
                    }
                    else
                    {
                        int64_t idx = (i-1)*yz3 + j*z3 + k*3 + 1;
                        if (idx < 0) continue;
                        indices[7] = shared_indices[idx];
                    }
                }
                if(edges & 0x100)
                {
                    if(i == 0 || j == 0)
                    {
                        indices[8] = vertices.size() / 3;
                        mc_add_vertex(x, y, z, z_dz, 2, v[0], v[4], _isoValue, &vertices);
                    }
                    else
                    {
                        int64_t idx = (i-1)*yz3 + (j-1)*z3 + k*3 + 2;
                        if (idx < 0) continue;
                        indices[8] = shared_indices[idx];
                    }
                }
                if(edges & 0x200)
                {
                    if(j == 0)
                    {
                        indices[9] = vertices.size() / 3;
                        mc_add_vertex(x_dx, y, z, z_dz, 2, v[1], v[5], _isoValue, &vertices);
                    }
                    else
                    {
                        int64_t idx = i*yz3 + (j-1)*z3 + k*3 + 2;
                        if (idx < 0) continue;
                        indices[9] = shared_indices[idx];
                    }
                }
                if(edges & 0x800)
                {
                    if(i == 0)
                    {
                        indices[11] = vertices.size() / 3;
                        mc_add_vertex(x, y_dy, z, z_dz, 2, v[3], v[7], _isoValue, &vertices);
                    }
                    else
                    {
                        int64_t idx = (i-1)*yz3 + j*z3 + k*3 + 2;
                        if (idx < 0) continue;
                        indices[11] = shared_indices[idx];
                    }
                }

                int tri;
                int* triangle_table_ptr = triangle_table[cubeindex];
                for(int m=0; tri = triangle_table_ptr[m], tri != -1; ++m)
                    polygons.push_back(indices[tri]);
            }
        }
    }


    // Convert the triangles to vertices
    for (size_t i = 0; i < vertices.size()/3; i++)
    {
        size_t x = vertices[i*3+0];
        size_t y = vertices[i*3+1];
        size_t z = vertices[i*3+2];

        Vector3f v(x, y, z);
        verticesA.push_back(v);
    }

    // Convert the triangles to vertices
    for (size_t i = 0; i < polygons.size() / 3; i++)
    {
        size_t i1 = int(polygons[i*3+0]);
        size_t i2 = int(polygons[i*3+1]);
        size_t i3 = int(polygons[i*3+2]);

        Triangle t;
        t[0] = i1; // verticesA[index + 0];
        t[1] = i2; // verticesA[index + 1];
        t[2] = i3; // verticesA[index + 2];

        trianglesA.push_back(t);

    }


     delete [] shared_indices;

    return  new Mesh(verticesA, trianglesA);
    return nullptr;
}

//static Vector<float> voxels(65*65*65);

//static void generate_voxels()
//{
//    Noise2D n2d(0);
//    for (int z = 0; z < 65; z++) {
//        for (int y = 0; y < 65; y++) {
//            for (int x = 0; x < 65; x++) {
//                const float fy = (float)y / 65.0f;
//                const int offset = offset_3d({x, y, z}, Vec3i(65));
//                const float v = n2d.get(x / 16.0f, z / 16.0f) * 0.25f;
//                voxels[offset] = fy - 0.25f - v;
//            }}}
//}



//struct Vertex {
//        Vec3f position;
//        Vec3f normal;
//};

//static Vector<float> voxels(65*65*65);
//static Vector<Vertex> vertices;

//static void triangle(int a, int b, int c)
//{
//    Vertex &va = vertices[a];
//    Vertex &vb = vertices[b];
//    Vertex &vc = vertices[c];
//    const Vec3f ab = va.position - vb.position;
//    const Vec3f cb = vc.position - vb.position;
//    const Vec3f n = cross(cb, ab);
//    va.normal += n;
//    vb.normal += n;
//    vc.normal += n;
//}

//static inline int offset_3d(const Vec3i &p, const Vec3i &size)
////{
////    return (p.z * size.y + p.y) * size.x + p.x;
////}


struct IdPoint
{
    size_t id;
    float x, y, z;
};

typedef std::map<size_t, IdPoint> PointIdMapping;

Mesh* reindex(PointIdMapping& vertexMapping, Triangles &triangles)
{
    size_t nextId = 0;

    for (auto& pair : vertexMapping)
    {
        pair.second.id = nextId++;
    }
    for (auto& triangle : triangles)
    {
        triangle[0] = vertexMapping[triangle[0]].id;
        triangle[1] = vertexMapping[triangle[1]].id;
        triangle[2] = vertexMapping[triangle[2]].id;
    }

    size_t vertexCount = vertexMapping.size();
    Vertices vertices;
    vertices.resize(vertexCount);
    size_t index = 0;
    for (const auto& pair : vertexMapping)
    {
        vertices[index][0] = pair.second.x;
        vertices[index][1] = pair.second.y;
        vertices[index][2] = pair.second.z;
        ++index;
    }

    size_t triangleCount = triangles.size();
    auto triangleIndices = new size_t[triangleCount * 3];
    index = 0;
    for (const auto& triangle : triangles)
    {
        triangleIndices[index * 3] = triangle[0];
        triangleIndices[index * 3 + 1] = triangle[1];
        triangleIndices[index * 3 + 2] = triangle[2];
        ++index;
    }

    Mesh* mesh = new Mesh(vertices, triangles);
    return mesh;
}


Mesh* MarchingCubes::generateMesh()
{
    std::vector<int> indices;

    // Build the mesh
    Vertices vertices;
    Triangles triangles;

    for (int z = -5; z < _volume->getDepth() + 5; z++)
    {
        LOOP_PROGRESS(z, _volume->getDepth());

        for (int y = -5; y < _volume->getHeight() + 5; y++)
        {
            for (int x = -5; x < _volume->getWidth() + 5; x++)
            {
                // Get the value of the voxel (around a cube)
                const uint8_t voxelConfig[8] =
                {
                    _volume->getConfirmedValue(x, y, z),
                    _volume->getConfirmedValue(x + 1, y, z),
                    _volume->getConfirmedValue(x, y + 1, z),
                    _volume->getConfirmedValue(x + 1, y + 1, z),

                    _volume->getConfirmedValue(x, y, z + 1),
                    _volume->getConfirmedValue(x + 1, y, z + 1),
                    _volume->getConfirmedValue(x, y + 1, z + 1),
                    _volume->getConfirmedValue(x + 1, y + 1, z + 1)
                };

//                const float voxelConfig[8] =
//                {
//                    voxels[offset_3d({x,   y,   z},   Vec3i(65))],
//                    voxels[offset_3d({x+1, y,   z},   Vec3i(65))],
//                    voxels[offset_3d({x,   y+1, z},   Vec3i(65))],
//                    voxels[offset_3d({x+1, y+1, z},   Vec3i(65))],
//                    voxels[offset_3d({x,   y,   z+1}, Vec3i(65))],
//                    voxels[offset_3d({x+1, y,   z+1}, Vec3i(65))],
//                    voxels[offset_3d({x,   y+1, z+1}, Vec3i(65))],
//                    voxels[offset_3d({x+1, y+1, z+1}, Vec3i(65))],
//                };


                size_t tableIndex = 0;
//                                if (volume[z*xyWidth + y*xWidth + x] < isoLevel)
//                                    tableIndex |= 1;
//                                if (volume[z*xyWidth + (y + 1)*xWidth + x] < isoLevel)
//                                    tableIndex |= 2;
//                                if (volume[z*xyWidth + (y + 1)*xWidth + (x + 1)] < isoLevel)
//                                    tableIndex |= 4;
//                                if (volume[z*xyWidth + y*xWidth + (x + 1)] < isoLevel)
//                                    tableIndex |= 8;
//                                if (volume[(z + 1)*xyWidth + y*xWidth + x] < isoLevel)
//                                    tableIndex |= 16;
//                                if (volume[(z + 1)*xyWidth + (y + 1)*xWidth + x] < isoLevel)
//                                    tableIndex |= 32;
//                                if (volume[(z + 1)*xyWidth + (y + 1)*xWidth + (x + 1)] < isoLevel)
//                                    tableIndex |= 64;
//                                if (volume[(z + 1)*xyWidth + y*xWidth + (x + 1)] < isoLevel)
//                                    tableIndex |= 128;

                const int config_n =
                        ((voxelConfig[0] < _isoValue) << 0) |
                        ((voxelConfig[1] < _isoValue) << 1) |
                        ((voxelConfig[2] < _isoValue) << 2) |
                        ((voxelConfig[3] < _isoValue) << 3) |
                        ((voxelConfig[4] < _isoValue) << 4) |
                        ((voxelConfig[5] < _isoValue) << 5) |
                        ((voxelConfig[6] < _isoValue) << 6) |
                        ((voxelConfig[7] < _isoValue) << 7);

                if (config_n == 0 || config_n == 255)
                    continue;


                int edge_indices[12];
                auto do_edge = [&](int n_edge, uint8_t va, uint8_t vb, int axis, const Vector3f &base)
                {
                    if ((va < _isoValue) == (vb < _isoValue))
                        return;

                    Vector3f v = base;
                    v[axis] += va / (va - vb);
                    edge_indices[n_edge] = vertices.size();
                    vertices.push_back(v);
                };

                do_edge(0,  voxelConfig[0], voxelConfig[1], 0, Vector3f(x, y,   z));
                do_edge(1,  voxelConfig[2], voxelConfig[3], 0, Vector3f(x, y+1, z));
                do_edge(2,  voxelConfig[4], voxelConfig[5], 0, Vector3f(x, y,   z+1));
                do_edge(3,  voxelConfig[6], voxelConfig[7], 0, Vector3f(x, y+1, z+1));

                do_edge(4,  voxelConfig[0], voxelConfig[2], 1, Vector3f(x,   y, z));
                do_edge(5,  voxelConfig[1], voxelConfig[3], 1, Vector3f(x+1, y, z));
                do_edge(6,  voxelConfig[4], voxelConfig[6], 1, Vector3f(x,   y, z+1));
                do_edge(7,  voxelConfig[5], voxelConfig[7], 1, Vector3f(x+1, y, z+1));

                do_edge(8,  voxelConfig[0], voxelConfig[4], 2, Vector3f(x,   y,   z));
                do_edge(9,  voxelConfig[1], voxelConfig[5], 2, Vector3f(x+1, y,   z));
                do_edge(10, voxelConfig[2], voxelConfig[6], 2, Vector3f(x,   y+1, z));
                do_edge(11, voxelConfig[3], voxelConfig[7], 2, Vector3f(x+1, y+1, z));

                const uint64_t config = marching_cube_tris[config_n];
                const int n_triangles = config & 0xF;
                const int n_indices = n_triangles * 3;
                const int index_base = indices.size();


                int offset = 4;
                for (int i = 0; i < n_indices; i++)
                {
                    const int edge = (config >> offset) & 0xF;
                    indices.push_back(edge_indices[edge]);
                    offset += 4;
                }

                for (int i = 0; i < n_triangles; i++) {
                    Triangle t;
                    t[0] = indices[index_base+i*3+0];
                    t[1] = indices[index_base+i*3+1];
                    t[2] = indices[index_base+i*3+2];

                    triangles.push_back(t);
                }
            }
        }
    }

    std::cout << triangles.size() << "\n";
    std::cout << vertices.size()<< "\n";

        Mesh* mesh = new Mesh(vertices, triangles);

//        // Statistics
//        LOG_STATUS_IMPORTANT("Mesh Reconstruction with DMC Stats.");
//        LOG_STATS(_dmcGenerationTime);

        return mesh;
}


Mesh* MarchingCubes::generateMeshFromVolume(Volume* volume)
{
    // Reconstruct a watertight mesh from the volume with DMC
    std::unique_ptr< MarchingCubes > workflowMC =
            std::make_unique< MarchingCubes >(volume);

    // Generate the DMC mesh
    return workflowMC->generateMesh();
}

}

//static inline int offset_3d_slab(const Vec3i &p, const Vec3i &size)
//{
//    return size.x * size.y * (p.z % 2) + p.y * size.x + p.x;
//}



//static void generate_geometry_smooth()
//{
//    static Vector<Vec3i> slab_inds(65*65*2);

//    for (int z = 0; z < 64; z++) {
//        for (int y = 0; y < 64; y++) {
//            for (int x = 0; x < 64; x++) {
//                const Vec3i p(x, y, z);
//                const float voxelConfig[8] = {
//                    voxels[offset_3d({x,   y,   z},   Vec3i(65))],
//                    voxels[offset_3d({x+1, y,   z},   Vec3i(65))],
//                    voxels[offset_3d({x,   y+1, z},   Vec3i(65))],
//                    voxels[offset_3d({x+1, y+1, z},   Vec3i(65))],
//                    voxels[offset_3d({x,   y,   z+1}, Vec3i(65))],
//                    voxels[offset_3d({x+1, y,   z+1}, Vec3i(65))],
//                    voxels[offset_3d({x,   y+1, z+1}, Vec3i(65))],
//                    voxels[offset_3d({x+1, y+1, z+1}, Vec3i(65))],
//                };

//                const int config_n =
//                        ((voxelConfig[0] < 0.0f) << 0) |
//                        ((voxelConfig[1] < 0.0f) << 1) |
//                        ((voxelConfig[2] < 0.0f) << 2) |
//                        ((voxelConfig[3] < 0.0f) << 3) |
//                        ((voxelConfig[4] < 0.0f) << 4) |
//                        ((voxelConfig[5] < 0.0f) << 5) |
//                        ((voxelConfig[6] < 0.0f) << 6) |
//                        ((voxelConfig[7] < 0.0f) << 7);

//                if (config_n == 0 || config_n == 255)
//                    continue;

//                auto do_edge = [&](int n_edge, float va, float vb, int axis, const Vec3i &p) {
//                    if ((va < 0.0) == (vb < 0.0))
//                        return;

//                    Vec3f v = ToVec3f(p);
//                    v[axis] += va / (va - vb);
//                    slab_inds[offset_3d_slab(p, Vec3i(65))][axis] = vertices.length();
//                    vertices.append({v, Vec3f(0)});
//                };

//                if (p.y == 0 && p.z == 0)
//                    do_edge(0,  voxelConfig[0], voxelConfig[1], 0, Vec3i(x, y,   z));
//                if (p.z == 0)
//                    do_edge(1,  voxelConfig[2], voxelConfig[3], 0, Vec3i(x, y+1, z));
//                if (p.y == 0)
//                    do_edge(2,  voxelConfig[4], voxelConfig[5], 0, Vec3i(x, y,   z+1));
//                do_edge(3,  voxelConfig[6], voxelConfig[7], 0, Vec3i(x, y+1, z+1));

//                if (p.x == 0 && p.z == 0)
//                    do_edge(4,  voxelConfig[0], voxelConfig[2], 1, Vec3i(x,   y, z));
//                if (p.z == 0)
//                    do_edge(5,  voxelConfig[1], voxelConfig[3], 1, Vec3i(x+1, y, z));
//                if (p.x == 0)
//                    do_edge(6,  voxelConfig[4], voxelConfig[6], 1, Vec3i(x,   y, z+1));
//                do_edge(7,  voxelConfig[5], voxelConfig[7], 1, Vec3i(x+1, y, z+1));

//                if (p.x == 0 && p.y == 0)
//                    do_edge(8,  voxelConfig[0], voxelConfig[4], 2, Vec3i(x,   y,   z));
//                if (p.y == 0)
//                    do_edge(9,  voxelConfig[1], voxelConfig[5], 2, Vec3i(x+1, y,   z));
//                if (p.x == 0)
//                    do_edge(10, voxelConfig[2], voxelConfig[6], 2, Vec3i(x,   y+1, z));
//                do_edge(11, voxelConfig[3], voxelConfig[7], 2, Vec3i(x+1, y+1, z));

//                int edge_indices[12];
//                edge_indices[0]  = slab_inds[offset_3d_slab({p.x, p.y,   p.z  }, Vec3i(65))].x;
//                edge_indices[1]  = slab_inds[offset_3d_slab({p.x, p.y+1, p.z  }, Vec3i(65))].x;
//                edge_indices[2]  = slab_inds[offset_3d_slab({p.x, p.y,   p.z+1}, Vec3i(65))].x;
//                edge_indices[3]  = slab_inds[offset_3d_slab({p.x, p.y+1, p.z+1}, Vec3i(65))].x;
//                edge_indices[4]  = slab_inds[offset_3d_slab({p.x,   p.y, p.z  }, Vec3i(65))].y;
//                edge_indices[5]  = slab_inds[offset_3d_slab({p.x+1, p.y, p.z  }, Vec3i(65))].y;
//                edge_indices[6]  = slab_inds[offset_3d_slab({p.x,   p.y, p.z+1}, Vec3i(65))].y;
//                edge_indices[7]  = slab_inds[offset_3d_slab({p.x+1, p.y, p.z+1}, Vec3i(65))].y;
//                edge_indices[8]  = slab_inds[offset_3d_slab({p.x,   p.y,   p.z}, Vec3i(65))].z;
//                edge_indices[9]  = slab_inds[offset_3d_slab({p.x+1, p.y,   p.z}, Vec3i(65))].z;
//                edge_indices[10] = slab_inds[offset_3d_slab({p.x,   p.y+1, p.z}, Vec3i(65))].z;
//                edge_indices[11] = slab_inds[offset_3d_slab({p.x+1, p.y+1, p.z}, Vec3i(65))].z;

//                const uint64_t config = marching_cube_tris[config_n];
//                const int n_triangles = config & 0xF;
//                const int n_indices = n_triangles * 3;
//                const int index_base = indices.length();

//                int offset = 4;
//                for (int i = 0; i < n_indices; i++) {
//                    const int edge = (config >> offset) & 0xF;
//                    indices.append(edge_indices[edge]);
//                    offset += 4;
//                }
//                for (int i = 0; i < n_triangles; i++) {
//                    triangle(
//                                indices[index_base+i*3+0],
//                            indices[index_base+i*3+1],
//                            indices[index_base+i*3+2]);
//                }
//            }}}
//    for (Vertex &v : vertices)
//        v.normal = normalize(v.normal);
//}

//static void quad(bool flip, int ia, int ib, int ic, int id)
//{
//    if (flip)
//        std::swap(ib, id);

//    Vertex &va = vertices[ia];
//    Vertex &vb = vertices[ib];
//    Vertex &vc = vertices[ic];
//    Vertex &vd = vertices[id];

//    const Vec3f ab = va.position - vb.position;
//    const Vec3f cb = vc.position - vb.position;
//    const Vec3f n1 = cross(cb, ab);
//    va.normal += n1;
//    vb.normal += n1;
//    vc.normal += n1;

//    const Vec3f ac = va.position - vc.position;
//    const Vec3f dc = vd.position - vc.position;
//    const Vec3f n2 = cross(dc, ac);
//    va.normal += n2;
//    vc.normal += n2;
//    vd.normal += n2;

//    indices.append(ia);
//    indices.append(ib);
//    indices.append(ic);

//    indices.append(ia);
//    indices.append(ic);
//    indices.append(id);
//}


//static void generate_geometry_naive_surface_nets()
//{
//    static Vector<int> inds(65*65*2);

//    for (int z = 0; z < 64; z++) {
//        for (int y = 0; y < 64; y++) {
//            for (int x = 0; x < 64; x++) {
//                const Vec3i p(x, y, z);
//                const float voxelConfig[8] = {
//                    voxels[offset_3d({x,   y,   z},   Vec3i(65))],
//                    voxels[offset_3d({x+1, y,   z},   Vec3i(65))],
//                    voxels[offset_3d({x,   y+1, z},   Vec3i(65))],
//                    voxels[offset_3d({x+1, y+1, z},   Vec3i(65))],
//                    voxels[offset_3d({x,   y,   z+1}, Vec3i(65))],
//                    voxels[offset_3d({x+1, y,   z+1}, Vec3i(65))],
//                    voxels[offset_3d({x,   y+1, z+1}, Vec3i(65))],
//                    voxels[offset_3d({x+1, y+1, z+1}, Vec3i(65))],
//                };

//                const int config_n =
//                        ((voxelConfig[0] < 0.0f) << 0) |
//                        ((voxelConfig[1] < 0.0f) << 1) |
//                        ((voxelConfig[2] < 0.0f) << 2) |
//                        ((voxelConfig[3] < 0.0f) << 3) |
//                        ((voxelConfig[4] < 0.0f) << 4) |
//                        ((voxelConfig[5] < 0.0f) << 5) |
//                        ((voxelConfig[6] < 0.0f) << 6) |
//                        ((voxelConfig[7] < 0.0f) << 7);

//                if (config_n == 0 || config_n == 255)
//                    continue;

//                Vec3f average(0);
//                int average_n = 0;
//                auto do_edge = [&](float va, float vb, int axis, const Vec3i &p) {
//                    if ((va < 0.0) == (vb < 0.0))
//                        return;

//                    Vec3f v = ToVec3f(p);
//                    v[axis] += va / (va - vb);
//                    average += v;
//                    average_n++;
//                };

//                do_edge(voxelConfig[0], voxelConfig[1], 0, Vec3i(x, y,     z));
//                do_edge(voxelConfig[2], voxelConfig[3], 0, Vec3i(x, y+1,   z));
//                do_edge(voxelConfig[4], voxelConfig[5], 0, Vec3i(x, y,     z+1));
//                do_edge(voxelConfig[6], voxelConfig[7], 0, Vec3i(x, y+1,   z+1));
//                do_edge(voxelConfig[0], voxelConfig[2], 1, Vec3i(x,   y,   z));
//                do_edge(voxelConfig[1], voxelConfig[3], 1, Vec3i(x+1, y,   z));
//                do_edge(voxelConfig[4], voxelConfig[6], 1, Vec3i(x,   y,   z+1));
//                do_edge(voxelConfig[5], voxelConfig[7], 1, Vec3i(x+1, y,   z+1));
//                do_edge(voxelConfig[0], voxelConfig[4], 2, Vec3i(x,   y,   z));
//                do_edge(voxelConfig[1], voxelConfig[5], 2, Vec3i(x+1, y,   z));
//                do_edge(voxelConfig[2], voxelConfig[6], 2, Vec3i(x,   y+1, z));
//                do_edge(voxelConfig[3], voxelConfig[7], 2, Vec3i(x+1, y+1, z));

//                const Vec3f v = average / Vec3f(average_n);
//                inds[offset_3d_slab(p, Vec3i(65))] = vertices.length();
//                vertices.append({v, Vec3f(0)});

//                const bool flip = voxelConfig[0] < 0.0f;
//                if (p.y > 0 && p.z > 0 && (voxelConfig[0] < 0.0f) != (voxelConfig[1] < 0.0f)) {
//                    quad(flip,
//                         inds[offset_3d_slab(Vec3i(p.x, p.y,   p.z),   Vec3i(65))],
//                            inds[offset_3d_slab(Vec3i(p.x, p.y,   p.z-1), Vec3i(65))],
//                            inds[offset_3d_slab(Vec3i(p.x, p.y-1, p.z-1), Vec3i(65))],
//                            inds[offset_3d_slab(Vec3i(p.x, p.y-1, p.z),   Vec3i(65))]
//                            );
//                }
//                if (p.x > 0 && p.z > 0 && (voxelConfig[0] < 0.0f) != (voxelConfig[2] < 0.0f)) {
//                    quad(flip,
//                         inds[offset_3d_slab(Vec3i(p.x,   p.y, p.z),   Vec3i(65))],
//                            inds[offset_3d_slab(Vec3i(p.x-1, p.y, p.z),   Vec3i(65))],
//                            inds[offset_3d_slab(Vec3i(p.x-1, p.y, p.z-1), Vec3i(65))],
//                            inds[offset_3d_slab(Vec3i(p.x,   p.y, p.z-1), Vec3i(65))]
//                            );
//                }
//                if (p.x > 0 && p.y > 0 && (voxelConfig[0] < 0.0f) != (voxelConfig[4] < 0.0f)) {
//                    quad(flip,
//                         inds[offset_3d_slab(Vec3i(p.x,   p.y,   p.z), Vec3i(65))],
//                            inds[offset_3d_slab(Vec3i(p.x,   p.y-1, p.z), Vec3i(65))],
//                            inds[offset_3d_slab(Vec3i(p.x-1, p.y-1, p.z), Vec3i(65))],
//                            inds[offset_3d_slab(Vec3i(p.x-1, p.y,   p.z), Vec3i(65))]
//                            );
//                }
//            }}}
//    for (Vertex &v : vertices)
//        v.normal = normalize(v.normal);
//}

