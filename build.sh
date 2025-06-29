#!/bin/bash

# 设置库路径变量
OPENMESH_ROOT="./thirdparty/OpenMesh-8.1"
OPENMESH_LIB="${OPENMESH_ROOT}/build/Build/lib"
OPENMESH_INCLUDE="${OPENMESH_ROOT}/src"
EIGEN_PATH="./thirdparty/eigen-3.4.0"

# 验证OpenMesh库文件存在
if [ ! -f "${OPENMESH_LIB}/libOpenMeshCore.so.8.1" ]; then
    echo "错误: 找不到 libOpenMeshCore.so.8.1"
    echo "请检查路径: ${OPENMESH_LIB}"
    exit 1
fi

# 验证Eigen头文件存在
if [ ! -d "${EIGEN_PATH}" ]; then
    echo "错误: 找不到 Eigen 头文件目录"
    echo "请检查路径: ${EIGEN_PATH}"
    exit 1
fi

# 设置运行时库路径
export LD_LIBRARY_PATH="${OPENMESH_LIB}:${LD_LIBRARY_PATH}"

# 编译函数：添加错误检查
compile_app() {
    local app_name=$1
    local source_file=$2
    
    echo "编译 ${app_name}..."
    g++ -std=c++11 ${source_file} LSD.cpp \
        -I . \
        -I "${OPENMESH_INCLUDE}" \
        -I "${EIGEN_PATH}" \
        -L "${OPENMESH_LIB}" \
        -lOpenMeshCore \
        -O2 \
        -o ${app_name} \
        -pthread \
        -Wl,-rpath,'$ORIGIN/thirdparty/OpenMesh-8.1/build/Build/lib'  # 设置rpath
    
    # 检查编译结果
    if [ $? -ne 0 ]; then
        echo "错误: ${app_name} 编译失败!"
        exit 1
    fi
    
    echo "${app_name} 编译成功!"
}

# 编译应用程序
compile_app "LSD-denoising_mt" "LSD-denoising_mt.cpp"
compile_app "LSD-Gdata_mt" "LSD-Gdata_mt.cpp"

# 添加其他需要编译的程序...
# compile_app "另一个程序" "源文件.cpp"

echo "所有程序编译完成!"
