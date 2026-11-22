# --------------------------------------------------------------------------- #
#    PBeFSCheckModule.cmake                                                   #
#                                                                             #
#    Checks whether the C++20 module pbefs [see src/PBeFS.cppm] can be        #
#    built, and if not turns PBEFS_BUILD_MODULE off with a warning that says  #
#    why, so that the build goes on with the headers alone, as when the       #
#    option is OFF. The module needs CMake 3.28 or later, the Ninja           #
#    generator (Ninja 1.11 or later) or a Visual Studio one, and Clang 16 or  #
#    later with clang-scan-deps or GCC 15 or later.                           #
#                                                                             #
#        include(PBeFSCheckModule)                                            #
#        pbefs_check_module()                                                 #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
# --------------------------------------------------------------------------- #

function(pbefs_check_module)
    if (NOT PBEFS_BUILD_MODULE)
        return()
    endif ()

    set(why "")
    if (CMAKE_VERSION VERSION_LESS 3.28)
        set(why "CMake ${CMAKE_VERSION} is older than 3.28")
    elseif (CMAKE_GENERATOR MATCHES "^Ninja")
        execute_process(COMMAND "${CMAKE_MAKE_PROGRAM}" --version
                OUTPUT_VARIABLE ninja OUTPUT_STRIP_TRAILING_WHITESPACE
                RESULT_VARIABLE res)
        if (NOT res EQUAL 0 OR ninja VERSION_LESS 1.11)
            set(why "Ninja ${ninja} is older than 1.11")
        endif ()
    elseif (NOT CMAKE_GENERATOR MATCHES "^Visual Studio")
        set(why "the generator ${CMAKE_GENERATOR} builds no modules (Ninja "
                "or Visual Studio is needed)")
    endif ()

    if (NOT why)
        set(id "${CMAKE_CXX_COMPILER_ID}")
        set(ver "${CMAKE_CXX_COMPILER_VERSION}")
        if (NOT ((id STREQUAL "Clang" AND ver VERSION_GREATER_EQUAL 16) OR
                 (id STREQUAL "GNU" AND ver VERSION_GREATER_EQUAL 15)))
            set(why "the compiler ${id} ${ver} does not build the module "
                    "(Clang 16 or GCC 15 or later is needed)")
        elseif (id STREQUAL "Clang" AND NOT CMAKE_CXX_COMPILER_CLANG_SCAN_DEPS)
            set(why "clang-scan-deps, which scans the dependencies of the "
                    "modules for Clang, is not found")
        endif ()
    endif ()

    if (why)
        string(CONCAT why ${why})
        message(WARNING "PBEFS_BUILD_MODULE is turned OFF: ${why}")
        set(PBEFS_BUILD_MODULE OFF PARENT_SCOPE)
    endif ()
endfunction()

# --------------------------------------------------------------------------- #
