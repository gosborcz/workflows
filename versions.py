import os
import json
import subprocess
from typing import Dict
from typing import Set


def find_command_aptget():
    line_found = os.popen(
        "zgrep 'Commandline: apt-get install -y' /var/log/apt/history.log | sed 's/--[^[:blank:]]*//g'").read()
    command_list = line_found.split("\n")[:-1]
    commands = []
    for cnt, command in enumerate(command_list):
        command_list[cnt] = command[32:]
        for el in command_list[cnt].split():
            commands.append(el)
    return commands


def find_command_apt():
    line_found = os.popen(
        "zgrep 'Commandline: apt install -y' /var/log/apt/history.log | sed 's/--[^[:blank:]]*//g'").read()
    command_list = line_found.split("\n")[:-1]
    commands = []
    for cnt, command in enumerate(command_list):
        command_list[cnt] = command[28:]
        for el in command_list[cnt].split():
            commands.append(el)
    return commands


def conc_apts(aptget_var, apt_var):
    names = []
    if aptget_var and apt_var:
        names = aptget_var + apt_var
    if not apt_var:
        names = aptget_var
    if not aptget_var:
        names = apt_var
    return names


def get_versions(names):
    versions = []
    for name in names:
        apt_cache_str = 'apt-cache policy ' + str(name).split('=')[0]
        apt_cache = os.popen(apt_cache_str).read()
        lines = apt_cache.split('\n')
        version = lines[4][5:-4]
        versions.append(version)

    return versions


def create_json(names, versions):
    apt_dict = {}
    for name, version in zip(names, versions):
        apt_dict[name] = version
    return apt_dict


def write_json(final_json, name):
    with open(name, 'w') as outfile:
        json.dump(final_json, outfile, indent=4)


def custom_tools_ver():
    if os.path.exists('/tools'):
        tools_list = os.listdir("/tools")
        versions_dict = {}
        for tool in tools_list:
            version = os.listdir("/tools/" + str(tool))[0]
            versions_dict[tool] = version
    else:
        versions_dict = {}
    return versions_dict


def tools_intelliseq_ver():
    if os.path.exists('/intelliseqtools'):
        tools_list = [f for f in os.listdir("/intelliseqtools") if f.endswith(('.py', '.sh', '.bash', 'r', 'R'))]
        versions_intelliseq_dict = {}
        if tools_list:
            for tool in tools_list:
                extension = tool.split('.')[-1]
                if extension == "py":
                    path_script = '/intelliseqtools/' + tool
                    cmd = ['python3', path_script, '--version']
                    version_out_line = str(subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0])
                    version = version_out_line.split()[-1][:-3]
                # condition for other extensions
                elif extension == "sh" or extension == "bash":
                    path_script = '/intelliseqtools/' + tool
                    cmd = ['bash', path_script, '--version']
                    version_out_line = str(subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0])
                    version = version_out_line.split()[-1][:-3]
                elif extension == "r" or extension == "R":
                    path_script = '/intelliseqtools/' + tool
                    cmd = ['Rscript', path_script, '--version']
                    version_out_line = str(subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0])
                    version = version_out_line.split()[-1][:-3]
                else:
                    version = ""
                versions_intelliseq_dict[tool] = version
    else:
        versions_intelliseq_dict = {}
    return versions_intelliseq_dict


def resource_ver():
    if os.path.exists('/resources'):
        resources = [d for d in os.listdir("/resources") if os.path.isdir(os.path.join("/resources", d))]
        resource_dict = {}
        for resource in resources:
            version = os.listdir("/resources/" + str(resource))[0]
            resource_dict[resource] = {}
            resource_dict[resource]["version"] = version
            if os.path.exists('/resources/{}/{}/subresources-versions/'.format(resource, version)):
                subresources = os.listdir('/resources/{}/{}/subresources-versions/'.format(resource, version))
                subresource_dict = {}
                for subresource in subresources:
                    sub_version = \
                        os.listdir(
                            '/resources/{}/{}/subresources-versions/'.format(resource, version) + str(subresource))[
                            0]
                    subresource_dict[subresource] = sub_version
                resource_dict[resource]["subresources"] = subresource_dict
            else:
                continue
    else:
        resource_dict = {}
    return resource_dict


def pip_installed_packages_versions(directory: str = '/root/.local/bin') -> Dict[str, str]:
    directory_content = os.listdir(directory)
    return {console_script:get_version_from_python_executable(console_script) for console_script in directory_content}


def get_version_from_python_executable(executable_name: str) -> str:
    cmd = [executable_name, '--version']
    p = subprocess.run(cmd, stdout=subprocess.PIPE, check=True)
    out_str = p.stdout.decode()
    return out_str.split()[-1]


def main():
    # aptget apt tools bco
    aptget_var = find_command_aptget()
    apt_var = find_command_apt()
    names = conc_apts(aptget_var, apt_var)
    versions = get_versions(names)
    apt_dict = create_json(names, versions)
    write_json(apt_dict, 'software_bco_apts.bco')

    # custom tools bco
    versions_dict = custom_tools_ver()
    write_json(versions_dict, 'software_bco_custom.bco')

    # intelliseq tools bco
    versions_intelliseq_dict = tools_intelliseq_ver()
    versions_intelliseq_dict.update(pip_installed_packages_versions())
    write_json(versions_intelliseq_dict, 'software_bco_intelliseq.bco')

    # resources bco

    resource_dict = resource_ver()
    write_json(resource_dict, 'datasource_bco_recources.bco')


if __name__ == '__main__': main()
