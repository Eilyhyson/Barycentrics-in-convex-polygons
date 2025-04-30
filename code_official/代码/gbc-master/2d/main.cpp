// Copyright Dmitry Anisimov danston@ymail.com (c) 2016-2017.

// Some examples on how to use the coordinate classes.

// README:
/*

    A quick way to incorporate my code in your code is to implement a wrapper that converts data
    stored in your data structures to my data structures. As a quick example:

    std::vector<VertexR2> square(4);

    square[0] = VertexR2(YourPointClass.x(), YourPointClass.y());
    square[1] = VertexR2(YourPointClass.x(), YourPointClass.y());
    square[2] = VertexR2(YourPointClass.x(), YourPointClass.y());
    square[3] = VertexR2(YourPointClass.x(), YourPointClass.y());

    VertexR2 centre(YourPointClass.x(), YourPointClass.y());

    std::vector<double> result;

    MyBarycentricClassR2 bary(square);
    bary.compute(centre, result);

    Below you can find other examples. For more information see README.md.

*/
#include <string>
#include <fstream>
#include <iostream>
#include <chrono>
// Local includes.
#include "./coords/MeanValueR2.hpp"
#include "./coords/HarmonicR2.hpp"
#include "./coords/PointwiseIterativeR2.h"
#include "./coords/NapoleonMvcR2.hpp"

#include "./extra/MeshR2.hpp"
#include "./extra/VertexR2.hpp"
#include "./extra/TriangulatorR2.hpp"
#include "./extra/AllCoordinatesR2.hpp"
#include "./extra/TestCoordinatesR2.hpp"
using namespace gbc;

bool readData(const std::string& path, std::vector<VertexR2> & poly) {
    std::ifstream ifs(path);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open file: " << path << std::endl;
        return false;
    }
    int rows, cols;
    ifs >> rows >> cols;
    // std::cout << rows << ' ' << cols << std::endl;
    // Resize the matrix to the appropriate size
    double x, y;
    for (int i = 0; i < rows; ++i) {
        ifs >> x >> y;
        poly.push_back(VertexR2(x, y));
    }
    ifs.close();
    return true;
}

bool readData(const std::string& path, std::vector<Face> & faces) {
    std::ifstream ifs(path);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open file: " << path << std::endl;
        return false;
    }
    int rows, cols;
    ifs >> rows >> cols;
    // Resize the matrix to the appropriate size
    double v1, v2, v3;
    for (int i = 0; i < rows; ++i) {
        ifs >> v1 >> v2 >> v3;
        faces.push_back(Face(v1, v2, v3));
    }
    ifs.close();
    return true;
}

bool readData(const std::string& path, std::vector<bool> & poly) {
    std::ifstream ifs(path);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open file: " << path << std::endl;
        return false;
    }
    int rows;
    ifs >> rows;
    // Resize the matrix to the appropriate size
    int x;
    for (int i = 0; i < rows; ++i) {
        ifs >> x;
        poly.push_back(x);
    }
    ifs.close();
    return true;
}

void load_imformation(const std::string path, std::vector<VertexR2> & poly, 
        std::vector<VertexR2> & x, std::vector<bool> & vertex_markers,
        std::vector<Face> & faces)      
{
    std::string path_poly = path + "polygon";
    // 读取polygon
    if (readData(path_poly, poly))
    {
        std::cout << "成功加载多边形" << std::endl;
    }

    std::string path_x = path + "x";
    // 读取x
    if (readData(path_x, x))
    {
        std::cout << "成功加载采样点" << std::endl;
    }

    std::string path_vertex_markers = path + "vertex_markers";
    // 读取vertex_markers
    if (readData(path_vertex_markers, vertex_markers))
    {
        std::cout << "成功加载vertex_markers" << std::endl;
    }

    std::string path_faces = path + "triangles";
    // 读取faces
    if (readData(path_faces, faces))
    {
        std::cout << "成功加载faces" << std::endl;
    }
} 

void load_imformation(const std::string path, std::vector<VertexR2> & poly)      
{
    std::string path_poly = path + "_polygon";
    // 读取polygon
    if (readData(path_poly, poly))
    {
        std::cout << "成功加载多边形" << std::endl;
    }
} 

bool save_BC(const std::string& path, std::vector<VertexR2> & x) {
    std::ofstream outputFile(path, std::ios::out);  // 创建文件并打开以写入数据
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open file: " << path << std::endl;
        return false;
    }
    // 写入数据
    for (auto& item : x) {
        std::vector<double> & bc = item.b(); 
        for (auto& bc_i : bc)
        {
            outputFile << bc_i << ' ';
        }
        outputFile << std::endl;
    }
    // 关闭文件
    outputFile.close();
    return true;
}

bool save_imformation(const std::string& path, std::vector<Face> & faces) {
    std::ofstream outputFile(path, std::ios::out);  // 创建文件并打开以写入数据
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open file: " << path << std::endl;
        return false;
    }
    // 写入数据
    for (auto& item : faces) {
        outputFile << item.v[0] << ' ' << item.v[1] << ' ' << item.v[2] << std::endl;
    }
    // 关闭文件
    outputFile.close();
    return true;
}

bool save_imformation(const std::string& path, std::vector<VertexR2> & x) {
    std::ofstream outputFile(path, std::ios::out);  // 创建文件并打开以写入数据
    if (!outputFile.is_open()) {
        std::cerr << "Unable to open file: " << path << std::endl;
        return false;
    }
    // 写入数据
    for (auto& item : x) {
        outputFile << item.x() << ' ' << item.y() << std::endl;
    }
    // 关闭文件
    outputFile.close();
    return true;
}

// Examples.
int main() {
    std::cout << "加载数据" << std::endl;
    const std::string path = "/home/hyson/公共/code_official/代码/data/鸭子";
    const std::string save_path = "/home/hyson/公共/code_official/代码/gbc-master/out/Cat_0.5/";
    std::vector<VertexR2> poly;
    load_imformation(path, poly);

    const double edgeLength = 0.5;  //修改参数调整采样点的密度
    // Refine the polygon to create regular mesh.
    std::vector<VertexR2> refined;
    const size_t n = poly.size();
    for (size_t i = 0; i < n; ++i) {
        refined.push_back(poly[i]);
        const size_t ip = (i + 1) % n;
        const size_t numS = ceil((poly[ip] - poly[i]).length() / edgeLength);
        for (size_t j = 1; j < numS; ++j) {
            VertexR2 vert = poly[i] + (double(j) / double(numS)) * (poly[ip] - poly[i]);
            refined.push_back(vert);
        }
    }

    // Create mesh. 生成网格
    TriangulatorR2 tri(refined, edgeLength, true);
    tri.setPlanarGraph(true);

    std::vector<VertexR2> queries;
    std::vector<Face> faces;

    tri.generate(queries, faces);
    MeshR2 mesh;
    mesh.initialize(queries, faces);
    // Clean mesh from the polygon vertices.
    std::vector<VertexR2> cleaned;
    for (size_t i = 0; i < mesh.numVertices(); ++i) {
        if (mesh.vertices()[i].type == INTERIOR || mesh.vertices()[i].type == FLAT)
            cleaned.push_back(mesh.vertices()[i]);
    }

    // 保存网格数据
    if (save_imformation(save_path + "x", queries))
    {
        std::cout << "保存采样点" << ' ' << queries.size() << std::endl;
    }
    if (save_imformation(save_path + "faces", faces))
    {
        std::cout << "保存faces" << std::endl;
    }

    //计算中值坐标
    auto last = std::chrono::high_resolution_clock::now();
    MeanValueR2 mvc = MeanValueR2(poly);
    mvc.bc(queries);
    auto cur = std::chrono::high_resolution_clock::now();
    // // 计算时间差
    std::chrono::duration<double> elapsed = cur - last;
    // 输出结果
    std::cout << "计算均值坐标时间: " << elapsed.count() << " 秒" << std::endl;
    if (save_BC(save_path + "mvc_lambda", queries))
    {
        std::cout << "中值坐标计算成功" << std::endl;
    }

    //计算正中值坐标
    last = std::chrono::high_resolution_clock::now();
    PositiveMeanValueR2 pmvc = PositiveMeanValueR2(poly);
    pmvc.bc(queries);
    cur = std::chrono::high_resolution_clock::now();
    // 计算时间差
    elapsed = cur - last;
    // 输出结果
    std::cout << "计算正中值坐标时间: " << elapsed.count() << " 秒" << std::endl;
    if (save_BC(save_path + "pmvc_lambda", queries))
    {
        std::cout << "正中值坐标计算成功" << std::endl;
    }

    //计算最大熵
    last = std::chrono::high_resolution_clock::now();
    MaximumEntropyR2 mec = MaximumEntropyR2(poly);
    mec.bc(queries);
    cur = std::chrono::high_resolution_clock::now();
    // 计算时间差
    elapsed = cur - last;
    // 输出结果
    std::cout << "最大熵坐标时间: " << elapsed.count() << " 秒" << std::endl;
    if (save_BC(save_path + "mec_lambda", queries))
    {
        std::cout << "最大熵坐标计算成功" << std::endl;
    }

     //计算调和坐标
    last = std::chrono::high_resolution_clock::now();
    HarmonicR2 hmc = HarmonicR2(poly);
    hmc.setMesh(queries, faces);
    hmc.bc(queries);
    cur = std::chrono::high_resolution_clock::now();
    // 计算时间差
    elapsed = cur - last;
    // 输出结果
    std::cout << "计算调和坐标时间: " << elapsed.count() << " 秒" << std::endl;
    if (save_BC(save_path + "hmc_lambda", queries))
    {
        std::cout << "调和坐标计算成功" << std::endl;
    }

    // 计算局部坐标
    last = std::chrono::high_resolution_clock::now();
    LocalR2 lbc = LocalR2(poly);
    lbc.setEdgeLength(edgeLength);
    lbc.setMesh(queries, faces);
    lbc.bc(queries);
    cur = std::chrono::high_resolution_clock::now();
    // 计算时间差
    elapsed = cur - last;
    // 输出结果
    std::cout << "计算局部坐标时间: " << elapsed.count() << " 秒" << std::endl;
    if (save_BC(save_path + "lbc_lambda", queries))
    {
        std::cout << "局部坐标计算成功" << std::endl;
    }


    // // 计算PIC
    // last = std::chrono::high_resolution_clock::now();
    // PointwiseIterativeR2 pic = PointwiseIterativeR2(poly);
    // pic.bc(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // // 计算时间差
    // elapsed = cur - last;
    // // 输出结果
    // std::cout << "计算pic坐标时间: " << elapsed.count() << " 秒" << std::endl;
    // if (save_BC(save_path + "pic_lambda", queries))
    // {
    //     std::cout << "pic坐标计算成功" << std::endl;
    // }
    // gbc::BarycentricPropertiesR2 test = BarycentricPropertiesR2(poly, 1e-8, true);
    // test.check_linearReproduction(queries);
    // test.check_partitionOfUnity(queries);

    //计算拿破仑三角形推到的mvc
    // last = std::chrono::high_resolution_clock::now();
    // NapoleonMvcR2 nmvc = NapoleonMvcR2(poly);
    // nmvc.bc(queries);
    // cur = std::chrono::high_resolution_clock::now();
    // // 计算时间差
    // elapsed = cur - last;
    // // 输出结果
    // std::cout << "计算nmvc坐标时间: " << elapsed.count() << " 秒" << std::endl;
    // if (save_BC(save_path + "nmvc_lambda", queries))
    // {
    //     std::cout << "nmvc坐标计算成功" << std::endl;
    // }
}

