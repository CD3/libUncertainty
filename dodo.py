import pathlib
import subprocess


def task_install_dependencies():
    return {
        "actions": ["conan install . --build missing"],
        "file_dep": ["conanfile.txt", __file__],
    }


def task_configure():
    return {
        "actions": ["cmake --preset conan-default"],
        "file_dep": ["CMakeLists.txt", __file__],
    }


def task_build():
    return {
        "actions": ["cmake --build --preset conan-release"],
        "file_dep": ["build/CMakeCache.txt", __file__],
    }


def task_test():
    def run():
        exes = pathlib.Path("build").rglob("libUncertainty_CatchTests")
        for exe in exes:
            subprocess.run(str(exe.relative_to("build")), shell=True, cwd="build")

    return {"actions": [run], "task_dep": ["build"]}
