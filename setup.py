#----------------------------------------------------------------------
from setuptools import find_packages, setup, Extension, Command
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install
# 
from distutils.command.build import build
#----------------------------------------------------------------------
import os
# OOP way to handle paths in python. Included in Standard Library from python 3.6
from pathlib import Path
# To count processors basically
import psutil
# to run shell commands
import subprocess
# to read json files
import json
#----------------------------------------------------------------------
# Custom modules
#from commander3.python import Downloader, Checker 
from pysetup import Downloader, Checker 
#----------------------------------------------------------------------

#print(python)
#Downloader.myprint()

'''
# Creating a custom command
class TestCommands(Command):
    """
    The "description" and "user_options" are mandatory!
    description is the description of the command as it appears after running python setup.py --help-commands
    user_options is the custom command-line arguments 
    """
    description = "My custom test command just to see if it will be visible with python setup.py --help-commands"
    user_options = [("myoption=", None, "Here is some description coming?")]
    #print(user_options)
    #print(type(user_options))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        print("You just run a custom command! Hooray :)")
'''
"""
Description: This class should override standard distutils (not even setup tools!) command.
Link to the source:
https://github.com/python/cpython/blob/master/Lib/distutils/command/build.py
Note: We do not extend the ditutils.build but instead creating a custom command. 
"""
class BuildCommand(Command):
    # Options for building the project. Format: (option name, short version of option name, description) 
    # Hephyns will be added by python automatically, i.e. parallel, means we use:
    # --parallel=4 or -j 4 to get a value of 4
    user_options = [
        ('parallel=', 'j', "Number of processors to use to compile the project."),
        ("build-dir=", None, "Name of the directory where to build the project."),
        ("download-dir=", None, "Name of the directory where to download (sub)projects."),
        ("install-dir=", None, "Name of the directory where to install (sub)projects.")
        ]
    # Command description
    description = "Overwritten distutils 'build' command."

    """
    Description: This method initialized command line options and sets their default values.
    The values are string values, so we convert them into strings here
    """
    def initialize_options(self):
        self.parallel     = psutil.cpu_count(logical=False) 
        self.build_dir    = Path("build").resolve() 
        self.download_dir = Path().joinpath(self.build_dir, "downloads").resolve()
        self.install_dir  = Path().joinpath(self.build_dir, "install").resolve()

    """
    Description: This method finilizes the process before running "build" command. 
    Aka, postprocessing options. 
    """
    def finalize_options(self):
        # Creating directories if they do not exist
        if not Path.is_dir(self.build_dir):
            Path.mkdir(self.build_dir)
        if not Path.is_dir(self.install_dir):
            Path.mkdir(self.install_dir)
        if not Path.is_dir(self.download_dir):
            Path.mkdir(self.download_dir)
        # Converting vars into string values as they should be
        self.build_dir    = str(self.build_dir) 
        self.download_dir = str(self.download_dir)
        self.install_dir  = str(self.install_dir)
        # Path to root of clone Commander3 repo
        #self.root_dir = Path(__file__).resolve().parent
        #self.project_dir = Path().joinpath(self.root_dir, self.build_dir)

        # Updating the environment variables safely:
        # Make a copy of the environment and pass is to the childen.
        self.comm3_env = os.environ.copy()
        self.comm3_env["PATH"] = f"{self.install_dir}/bin:{os.environ.get('PATH')}"
        self.comm3_env["LD_LIBRARY_PATH"] = f"{self.install_dir}/lib:{os.environ.get('LD_LIBRARY_PATH')}"
        # CMake commands
        self.cmake_configure_commands = ["cmake", 
                f"-DCMAKE_INSTALL_PREFIX={self.install_dir}",
                ".."
                ]

    """
    Description: This method actually runs the "build" command
    """
    def run(self):
        """
        print(f"-- ============================================")
        print(f"-- Start analyzing. please, wait a moment")
        result = subprocess.run(self.cmake_configure_commands,
                shell=False, check=True, env=self.comm3_env, cwd=self.build_dir)#,
        #        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        """
        print("THERE IS NO BUILD COMMAND AS OF NOW!")
        print("PLEASE, RUN `python setup.py install` to install Commander3")
        print("EXITING...")
        exit()
        #print("Running custom build command")
        #print(self.parallel)
        #print(self.build_dir)
        #print(self.download_dir)
        #print(self.install_dir)
        #build.run(self)

# TODO: probably add everything from BuildCommand to InstallCommand to have one-line command
class Installer(Command):
    # Options for building the project. Format: (option name, short version of option name, description) 
    # Hephyns will be added by python automatically, i.e. parallel, means we use:
    # --parallel=4 or -j 4 to get a value of 4
    user_options = [
        ('parallel=', 'j', "Number of processors to use to compile the project."),
        ("build-dir=", None, "Name of the directory where to build the project."),
        ("download-dir=", None, "Name of the directory where to download (sub)projects."),
        ("install-dir=", None, "Name of the directory where to install (sub)projects.")
        ]
    # Command description
    description = "Overwritten distutils 'build' command."

    """
    Description: This method initialized command line options and sets their default values.
    The values are string values, so we convert them into strings here
    """
    def initialize_options(self):
        self.parallel     = psutil.cpu_count(logical=False) 
        self.build_dir    = Path("build").resolve() 
        self.download_dir = Path().joinpath(self.build_dir, "downloads").resolve()
        self.install_dir  = Path().joinpath(self.build_dir, "install").resolve()

    """
    Description: This method finilizes the process before running "build" command. 
    Aka, postprocessing options. 
    """
    def finalize_options(self):
        # Creating directories if they do not exist
        if not Path.is_dir(self.build_dir):
            Path.mkdir(self.build_dir)
        if not Path.is_dir(self.install_dir):
            Path.mkdir(self.install_dir)
        if not Path.is_dir(self.download_dir):
            Path.mkdir(self.download_dir)
        # Path to root of clone Commander3 repo
        #self.root_dir = Path(__file__).resolve().parent
        #self.project_dir = Path().joinpath(self.root_dir, self.build_dir)

        # Updating the environment variables safely:
        # Make a copy of the environment and pass is to the childen.
        self.comm3_env = os.environ.copy()
        self.comm3_env["PATH"] = f"{self.install_dir}/bin:{os.environ.get('PATH')}"
        self.comm3_env["LD_LIBRARY_PATH"] = f"{self.install_dir}/lib:{os.environ.get('LD_LIBRARY_PATH')}"

    """
    Description: This method actually runs the "build" command
    """
    def run(self):
        print(f"-- ============================================")
        print(f"-- Start analyzing. please, wait a moment")
        #----------------------------------------------------------------------
        # Getting prerequisites
        root_dir = Path(__file__).resolve().parent
        prerequisites = str(root_dir.joinpath("prerequisites.json"))
        # read file
        with open(prerequisites, 'r') as myfile:
            data = myfile.read()
        # parse file
        prerequisites = json.loads(data)
        #----------------------------------------------------------------------
        # Getting GNU signatures for verification
        # Downloading GNU GPG keys
        gnu_signatures_source_url = "https://ftp.gnu.org/gnu/gnu-keyring.gpg"
        gnu_signatures_file = "gnu-keyring.gpg"
        full_path_to_gnu_signatures_file = self.download_dir / gnu_signatures_file
        Downloader.perform_download(gnu_signatures_source_url, 
                gnu_signatures_file, full_path_to_gnu_signatures_file)
        #----------------------------------------------------------------------
        # We get compilers with which we will compile commander
        num_procs = int(self.parallel)
        compilers_to_use, mpi_to_use = Checker.check_prerequisites(prerequisites, 
                self.comm3_env, self.download_dir, self.install_dir, num_procs, 
                full_path_to_gnu_signatures_file)
        print(f"--------------------------------------------")
        print(f"-- The following compilers will be used {compilers_to_use}")
        print(f"-- The following mpi wrappers will be used {mpi_to_use}")
        print(f"--------------------------------------------")
        #----------------------------------------------------------------------
        # Converting vars into string values as they should be
        #self.build_dir    = str(self.build_dir) 
        #self.download_dir = str(self.download_dir)
        #self.install_dir  = str(self.install_dir)
        #exit()
        #----------------------------------------------------------------------
        # CMake commands
        cmake_configure_command = ["cmake", 
                f"-DCMAKE_INSTALL_PREFIX={self.install_dir}",
                f"-DCMAKE_C_COMPILER={compilers_to_use[0]}",
                f"-DCMAKE_CXX_COMPILER={compilers_to_use[1]}",
                f"-DCMAKE_Fortran_COMPILER={compilers_to_use[2]}",
                f"-DMPI_C_COMPILER={mpi_to_use[0]}",
                f"-DMPI_CXX_COMPILER={mpi_to_use[1]}",
                f"-DMPI_Fortran_COMPILER={mpi_to_use[2]}",
                ".."
                ]
        cmake_configure_project = subprocess.run(cmake_configure_command,
                shell=False, check=True, env=self.comm3_env, cwd=self.build_dir)#,
        #        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        
        cmake_install_command = ["cmake", "--build", ".", "--target", "install", "-j", f"{num_procs}"]
        cmake_build_project = subprocess.run(cmake_install_command,
                shell=False, check=True, env=self.comm3_env, cwd=self.build_dir)#,
        #        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        #----------------------------------------------------------------------
        #print("Running custom build command")
        #print(self.parallel)
        #print(self.build_dir)
        #print(self.download_dir)
        #print(self.install_dir)
        #build.run(self)

#----------------------------------------------------------------------
"""
Description: Override the base extension class for obstructing setuptools from 
building the project for us. All subprojects will be compiled with either python
from source using custom scripts or CMake.
"""
class CMakeExtension(Extension):
    def __init__(self, name, sources=[]):
        super().__init__(name = name, sources = sources)

"""
Description: Building package using CMake
"""
class CMakeBuild(build_ext):

    def run(self):
        super().run()
#----------------------------------------------------------------------
# Parsing README.md file
readme_file = Path(__file__).parent.absolute()/"README.md"
def get_project_description(readme_file):

    pass
#----------------------------------------------------------------------
"""
Description: Override the base extension class for obstructing setuptools from 
building the project for us. All subprojects will be compiled with either python
from source using custom scripts or CMake.
"""
class CMakeExtension(Extension):
    def __init__(self, name, sources=[]):
        super().__init__(name = name, sources = sources)

"""
Description: Building package using CMake
"""
class CMakeBuild(build_ext):

    def run(self):
        super().run()
#----------------------------------------------------------------------
# Parsing README.md file
readme_file = Path(__file__).parent.absolute()/"README.md"
def get_project_description(readme_file):
    with open(readme_file, "r", encoding="utf-8") as readme_file:
        project_description = readme_file.read()
#----------------------------------------------------------------------
# Creating configure options for setuptools to install the project
setup_config_args = {}
# Populating setuptools options
setup_config_args['name'] = "commander3-cmb" 
setup_config_args['description'] = "Optimal Monte-carlo Markov chAiN Driven EstimatoR" 
setup_config_args['long_description'] = get_project_description(readme_file)
setup_config_args['long_description_content_type'] = "text/markdown"
setup_config_args['author'] = "Maksym Brilenkov" 
setup_config_args['author_email'] = "maksym.brilenkov@astro.uio.no"
setup_config_args['license'] = "GPLv3" 
setup_config_args['url'] = "https://github.com/Cosmoglobe/Commander" 
# TODO: make a unified version file, e.g. versions.cmake and parse it here as well
# This is just a placeholder now
setup_config_args['version'] = "1.0.0" 
setup_config_args['python_requires'] = ">=3.6.0"
setup_config_args['install_requires'] = [
        "cmake >= 3.17",
        "requests",
        # NOT gnupg or pretty_bad_protocols
        "python-gnupg",
        #"numpy",
        "packaging",
        "psutil"
        ] 
# TODO: it seems I need to overwrite an existing "build" and "install" commands to make them as CMake ones
# instead of build_ext
setup_config_args["cmdclass"] = {
        "build": BuildCommand,
        #"install": InstallCommand, 
        "install": Installer, 
        "build_ext": CMakeBuild 
        #"mycommand": TestCommands
        }
setup_config_args['classifiers'] = [
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Fortran",
        "Programming Language :: C",
        "Programming Language :: C++"
        ] 
#----------------------------------------------------------------------
# Initiating package setup
setup(**setup_config_args)
