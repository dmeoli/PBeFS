package("pbefs")
    set_homepage("https://github.com/dmeoli/PBeFS")
    set_description("Parallel best-first search in C++20: A* and branch-and-bound")
    set_license("LGPL-3.0")

    add_urls("https://github.com/dmeoli/PBeFS/releases/download/v$(version)/PBeFS-$(version).tar.gz",
             "https://github.com/dmeoli/PBeFS.git")
    add_versions("0.1.0", "0000000000000000000000000000000000000000000000000000000000000000")

    add_deps("cmake")
    add_deps("ads")

    if is_plat("linux", "bsd") then
        add_syslinks("pthread")
    end

    on_install(function (package)
        -- FastFlow, which would be downloaded at configure time, and MPI
        -- are left out
        local configs = {"-DPBEFS_WITH_FASTFLOW=OFF",
                         "-DPBEFS_WITH_MPI=OFF",
                         "-DPBEFS_BUILD_TESTS=OFF",
                         "-DPBEFS_BUILD_EXAMPLES=OFF",
                         "-DPBEFS_BUILD_TOOLS=OFF"}
        table.insert(configs, "-DCMAKE_BUILD_TYPE=" ..
                     (package:is_debug() and "Debug" or "Release"))
        table.insert(configs, "-DBUILD_SHARED_LIBS=" ..
                     (package:config("shared") and "ON" or "OFF"))
        import("package.tools.cmake").install(package, configs)
    end)

    on_test(function (package)
        assert(package:check_cxxsnippets({test = [[
            using namespace PBeFS_di_unipi_it;
            void test() {
                const EightPuzzle d(EightPuzzle::parse("7 2 4 5 0 6 8 3 1"));
                auto s = astar(d);
            }
        ]]}, {configs = {languages = "c++20"},
              includes = "PBeFS/PBeFS.h"}))
    end)
