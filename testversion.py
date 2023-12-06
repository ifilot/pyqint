# -*- coding: utf-8 -*-

import os
import re

ROOT = os.path.dirname(__file__)

def main():
    version_versionpy = get_version_versionpy()
    version_metayaml = get_version_metayaml()
    
    print('Version strings found:')
    print(version_versionpy)
    print(version_metayaml)

    try:
        for i in range(0,3):
             assert version_versionpy[i] == version_metayaml[i]
    except Exception as e:
        print(e)
        raise Exception('Invalid version strings encountered')

def get_version_projecttoml():
    """
    Extract the version string from the pyproject.toml file
    """
    pattern = re.compile(r'^version\s*=\s*"(\d+\.\d+.\d+)"\s*$')

    f = open(os.path.join(ROOT, 'pyproject.toml'))
    lines = f.readlines()
    for line in lines:
        match = re.match(pattern, line)
        if match:
            version = match.groups(1)[0]
            return [int(i) for i in version.split('.')]

    return None

def get_version_versionpy():
    """
    Extract the version string from the _version.py file
    """
    pattern = re.compile(r'^__version__\s*=\s*[\'"](\d+\.\d+.\d+)[\'"]\s*$')

    f = open(os.path.join(ROOT, 'pyqint', '_version.py'))
    lines = f.readlines()
    for line in lines:
        match = re.match(pattern, line)
        if match:
            version = match.groups(1)[0]
            return [int(i) for i in version.split('.')]

    return None

def get_version_metayaml():
    """
    Extract the version string from the meta.yaml file
    """
    pattern = re.compile(r'^\s*version\s*:\s*"(\d+\.\d+.\d+)"\s*$')

    f = open(os.path.join(ROOT, 'meta.yaml'))
    lines = f.readlines()
    for line in lines:
        match = re.match(pattern, line)
        if match:
            version = match.groups(1)[0]
            return [int(i) for i in version.split('.')]

    return None

if __name__ == '__main__':
    main()