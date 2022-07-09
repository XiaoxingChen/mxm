#!/usr/bin/python3
# General Build Script Template
import sys
import os
import shutil
import argparse
from time import time
sys.path.append( os.path.join(os.path.abspath(os.path.dirname(__file__)), 'scripts'))

UNIX_CMAKE_STR = 'cmake {} -B{} -G "Unix Makefiles"'
MSVC_CMAKE_STR = 'cmake {} -B{} -G "Visual Studio 16 2019"'
TEST_EXCUTABLE_NAME = 'test_main'

def handleAutoComplete():
    if sys.platform == 'linux':
        complete_cmd = 'complete -F _longopt {}'.format(os.path.basename(__file__))
        bashrc_path = os.path.expanduser('~/.bashrc')
        with open(bashrc_path) as f:
            if not complete_cmd in f.read():
                os.system('echo "{}" >> {}'.format(complete_cmd, bashrc_path))
    else:
        pass

class Dir():
    script_folder = os.path.abspath(os.path.dirname(__file__))
    build_root = os.path.join(script_folder, 'build')
    cmake_project = script_folder
    plot_script = os.path.join(script_folder, 'scripts', 'pixel_coordinate.py')
    lfs_yaml = os.path.join(script_folder, 'lfs.yaml')
    lfs_asset = os.path.join(build_root, 'assets')

class BuildTarget():
    def __init__(self):
        self.require_test = False
        self.require_build = False
        self.host_exec_relative = None
        self.exec_args = None
        self.cmake_options = ''

    @property
    def host_exec(self):
        if self.host_exec_relative is None:
            return None
        return os.path.join(self.build_dir, self.host_exec_relative)

    def updateFromArgs(self, args):
        self.require_test = args.test_all or args.__getattribute__(self.parser_test_var)
        platform_exe = args.__getattribute__(self.parser_run_var)

        self.require_build = True
        if self.require_test:
            self.host_exec_relative = os.path.join('tests', TEST_EXCUTABLE_NAME)
        elif platform_exe:
            self.host_exec_relative = platform_exe[0]
            self.exec_args = ' '.join(platform_exe[1:])
        elif args.__getattribute__(self.platform.replace('-', '_')):
            pass
        else:
            self.require_build = False

        if args.cmake_options:
            self.cmake_options = args.cmake_options

    def cmakeBaseOptions(self):
        CMAKE_OPTION_STR = ''
        if self.require_test:
            CMAKE_OPTION_STR += ' -DBUILD_TEST=1'
        if self.cmake_options:
            CMAKE_OPTION_STR += self.cmake_options
        return CMAKE_OPTION_STR

    def updateArgParser(self, parser):
        self.parser_build_cmd = ('--' + self.platform)

        self.parser_run_cmd = '--run' + ('' if self.platform == 'native' else ('-' + self.platform))
        self.parser_run_var = (self.platform) + '_exe'
        # self.parser_run_var.replace('-','_')

        self.parser_test_cmd = '--test' + ('' if self.platform == 'native' else ('-' + self.platform))
        self.parser_test_var = self.parser_test_cmd[2:].replace('-','_')

        parser.add_argument(self.parser_build_cmd, action='store_true', \
            help='Build ' + self.platform + ' platform')
        parser.add_argument(self.parser_run_cmd, dest=self.parser_run_var, \
            nargs='+', help='run ' + self.platform + ' executable')
        parser.add_argument(self.parser_test_cmd, action='store_true', \
            help='Run test on ' + self.platform + ' platform')


    def checkExecutableValid(self):
        if self.host_exec is None:
            return False
        if not os.path.isfile(self.host_exec):
            print('Error: {} executable {} not found!'.format(self.platform, self.host_exec))
            return False
        return True

    @property
    def build_dir(self):
        return os.path.join(Dir.build_root, self.platform)

class LinuxNativeTarget(BuildTarget):
    platform = 'native'
    def runBuild(self):
        if not self.require_build:
            return None
        cmake_generate_cmd = UNIX_CMAKE_STR.format(Dir.cmake_project, self.build_dir)
        cmake_generate_cmd += self.cmakeBaseOptions()

        os.system(cmake_generate_cmd)
        exit_code = os.system('cmake --build {} -- -j8'.format(self.build_dir))
        if exit_code != 0:
            exit(1)

    def runExecutable(self):
        if not self.checkExecutableValid():
            return None
        exit_code = os.system("{} {}".format(self.host_exec, self.exec_args))
        if exit_code != 0:
            exit(1)

class WindowsNativeTarget(BuildTarget):
    platform = 'native'
    def runBuild(self):
        if not self.require_build:
            return None
        cmake_generate_cmd = MSVC_CMAKE_STR.format(Dir.cmake_project, self.build_dir)
        cmake_generate_cmd += self.cmakeBaseOptions()

        os.system(cmake_generate_cmd)
        exit_code = os.system('cmake --build {}'.format(self.build_dir))
        if exit_code != 0:
            exit(1)

    def runExecutable(self):
        if not self.checkExecutableValid():
            return None
        exit_code = os.system("{} {}".format(self.host_exec, self.exec_args))
        if exit_code != 0:
            exit(1)

class AndroidTarget(BuildTarget):

    def runExecutable(self):
        if not self.checkExecutableValid():
            return None

        device_work_dir = '/data/local/tmp/home/chenxx/'
        device_exec = device_work_dir + self.host_exec_relative

        os.system('adb shell "mkdir -p {}"'.format(device_work_dir))
        os.system('adb push {} {}'.format(self.host_exec, device_exec))
        os.system('adb shell chmod +x {}'.format(device_exec))
        exit_code = os.system('adb shell LD_LIBRARY_PATH={} {} {}'.format(device_work_dir, device_exec, self.exec_args))
        if exit_code != 0:
            exit(1)

    def runBuild(self):
        if not self.require_build:
            return None

        tool_chain_file = os.path.join(os.environ['ANDROID_NDK'],'build','cmake','android.toolchain.cmake')
        CMAKE_ANDROID_STR = ' -DCMAKE_TOOLCHAIN_FILE={} -DANDROID_PLATFORM=android-26'.format(tool_chain_file)

        cmake_generate_cmd = UNIX_CMAKE_STR.format(Dir.cmake_project, self.build_dir)
        cmake_generate_cmd += CMAKE_ANDROID_STR
        if self.platform == 'android':
            cmake_generate_cmd += ' -DANDROID_ABI=arm64-v8a'
        elif self.platform == 'android-arm32':
            cmake_generate_cmd += ' -DANDROID_ABI=armeabi-v7a'
        cmake_generate_cmd += self.cmakeBaseOptions()

        os.system(cmake_generate_cmd)
        exit_code = os.system('cmake --build {} -- -j8'.format(self.build_dir))
        if exit_code != 0:
            exit(1)

class AndroidArm64Target(AndroidTarget):
    platform = 'android'

class AndroidArm32Target(AndroidTarget):
    platform = 'android-arm32'

def createParser():
    # initialize argparse options
    parser = argparse.ArgumentParser()

    parser.add_argument('--clean', action='store_true', help='Clean build folder')
    parser.add_argument('--plot', nargs='+', help='plot target image')
    parser.add_argument('--all', action='store_true', help='Compile for all platforms')
    parser.add_argument('--test-all', action='store_true', help='Run test on all platforms')
    parser.add_argument('--cmake-options', dest='cmake_options', help='additional cmake options. Put in quote, start with space: " -DCMAKE_BUILD_TYPE=Debug"')
    parser.add_argument('--sync-lfs', action='store_true', help='Synchronize large file storage')

    return parser


def run(build_script_folder=os.path.abspath(os.path.dirname(__file__))):
    parser = createParser()
    handleAutoComplete()

    native_target = WindowsNativeTarget() if sys.platform == 'win32' else LinuxNativeTarget()
    targets = [native_target, AndroidArm64Target(), AndroidArm32Target()]
    for target in targets:
        target.updateArgParser(parser)

    args = parser.parse_args()

    if args.clean:
        shutil.rmtree(Dir.build_root, ignore_errors=True)
        quit()

    if args.sync_lfs:
        WebDrive.sync(Dir.lfs_yaml, Dir.lfs_asset)

    t_start = time()

    build_target_cnt = 0
    for target in targets:
        target.updateFromArgs(args)
        target.runBuild()
        build_target_cnt += int(target.require_build)

    if 0 == build_target_cnt:
        targets[0].require_build = True
        targets[0].runBuild()

    t_finish = time()
    print("Total build time: {:.3f}s".format(t_finish - t_start))

    for target in targets:
        target.runExecutable()

    if args.plot:
        os.system('python3 {} {}'.format(Dir.plot_script, args.plot[0]))


if __name__ == "__main__":
    run()