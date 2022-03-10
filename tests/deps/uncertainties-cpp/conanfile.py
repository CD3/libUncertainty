from conans import ConanFile, tools
import os
from glob import glob


class UncertaintiesCppConan(ConanFile):
    name = "uncertainties-cpp"
    version = "0.0.1"
    topics = ("eigen", "algebra", "linear-algebra", "vector", "numerical")

    @property
    def _source_subfolder(self):
        return "source_subfolder"

    def source(self):
        source_url = "https://github.com/Gattocrucco/uncertainties-cpp"
        self.run(f"git clone {source_url} {self.name}")


    def package(self):
        self.copy("*.hpp",src="uncertainties-cpp",dst="include")

    def package_id(self):
        self.info.header_only()

    def package_info(self):
        self.cpp_info.name = "uncertainties-cpp"

