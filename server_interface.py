# Here are the functions that download, save, and
# analyze the data from completed dft calculations


def add_escapes(calculation_directory):
    for char in "( )":
        calculation_directory = calculation_directory.replace(char, "\\" + char)
    return calculation_directory


def upsync(
    calculation_directory_list, host_name, remote_directory, options="r", version="3"
):
    import os

    """Upload a set of files to a remote host.
    Parameters:
        calculation_directory_list : string
            The list of absolute paths to be uploaded
        host_name : string
            Name of the host computer (i.e. user@cypress1.tulane.edu)
        remote_directory : string
            The location on the remote file system to upload to.
        options : string
            rsync options (v - verbosity, r - recursive, etc.)
            defaults to r.
        version : int
            Choose your current local rsync version.
    """
    if len(calculation_directory_list) == 0:
        return
    # Clean paths
    calculation_directory_list = list(map(add_escapes, calculation_directory_list))
    remote_directory = add_escapes(remote_directory)
    sourceStr = " ".join(calculation_directory_list)
    destStr = "{}:{}".format(host_name, remote_directory)
    if "v" in options:
        print("rsync -{} {} {}".format(options, sourceStr, destStr))
    os.system("rsync -{} {} {}".format(options, sourceStr, destStr))


def downsync(calculation_directory_list, host_name, localDir, options="r", version=2):
    import os

    """Download a set of files from a remote host.
    Parameters:
        caclDirList: string
            The list of absolute paths to be downloaded
        host_name : string
            Name of the host computer (i.e. user@cypress1.tulane.edu)
        remote_directory : string
            The location on the local file system to download to.
        options : string
            rsync options (v - verbosity, r - recursive, etc)
            defaults to r.
        version : int | 2 or 3
            Choose your current local rsync version.
    """
    if len(calculation_directory_list) == 0:
        return
    # Clean paths
    calculation_directory_list = list(map(add_escapes, calculation_directory_list))
    localDir = add_escapes(localDir)
    if version == 2:
        sourceStr = host_name + ':"{}"'.format(" ".join(calculation_directory_list))
    else:
        sourceStr = host_name + ":{}".format(" :".join(calculation_directory_list))
    if "v" in options:
        print("rsync -{} {} {}".format(options, sourceStr, localDir))
    os.system("rsync -{} {} {}".format(options, sourceStr, localDir))
