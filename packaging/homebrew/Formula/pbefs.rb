class Pbefs < Formula
  desc "Parallel best-first search in C++20: A* and branch-and-bound"
  homepage "https://github.com/dmeoli/PBeFS"
  url "https://github.com/dmeoli/PBeFS/releases/download/v0.1.0/PBeFS-0.1.0.tar.gz"
  sha256 "0000000000000000000000000000000000000000000000000000000000000000"
  license "LGPL-3.0-only"
  head "https://github.com/dmeoli/PBeFS.git", branch: "main"

  depends_on "cmake" => [:build, :test]

  def install
    args = %w[
      -DBUILD_SHARED_LIBS=ON
      -DPBEFS_WITH_FASTFLOW=OFF
      -DPBEFS_WITH_MPI=OFF
      -DPBEFS_BUILD_TESTS=OFF
      -DPBEFS_BUILD_EXAMPLES=OFF
      -DPBEFS_BUILD_TOOLS=OFF
    ]
    system "cmake", "-S", ".", "-B", "build", *args, *std_cmake_args
    system "cmake", "--build", "build"
    system "cmake", "--install", "build"
  end

  test do
    (testpath/"CMakeLists.txt").write <<~CMAKE
      cmake_minimum_required(VERSION 3.21)
      project(use LANGUAGES CXX)
      find_package(PBeFS REQUIRED)
      add_executable(use use.cpp)
      target_link_libraries(use PRIVATE PBeFS::PBeFS)
    CMAKE
    (testpath/"use.cpp").write <<~CPP
      #include "PBeFS/PBeFS.h"
      using namespace PBeFS_di_unipi_it;
      int main() {
        const EightPuzzle d(EightPuzzle::parse("7 2 4 5 0 6 8 3 1"));
        return astar(d).cost == 26 ? 0 : 1;
      }
    CPP
    system "cmake", "-S", ".", "-B", "build", *std_cmake_args
    system "cmake", "--build", "build"
    system "./build/use"
  end
end
