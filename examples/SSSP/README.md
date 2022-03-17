# SSSP file naming conventions and standards

## SSSP configurations
The SSSP comes in various configurations, where a configuration is the combination of the parameters version, functional and protocol.
The currently valid values for each of these parameters are:

 * version:  1.0, 1.1, 1.1.1, 1.1.2
 * functional: PBE, PBEsol
 * protocol: efficiency, precision

Note that not every possible permutation of these properties necessarily corresponds to an existing SSSP configuration.
For example, there is no PBEsol configuration for SSSP 1.0.

Note also that for each minor version, only the files of the latest patch version are provided here.
The files of the previous patch versions (i.e. 1.1, 1.1.1) can be retrieved from the previous versions of this record. 

In order to use the SSSP configurations in AiiDA, it is recommended to install them through the [aiida-pseudo](https://aiida-pseudo.readthedocs.io/en/latest/) package.

## SSSP file naming conventions
Each SSSP configuration comes with two files, where the filename is formatted as `SSSP_{version}_{functional}_{protocol}.{extension}`:

 * A gzipped tarball (`.tar.gz`) containing the pseudo potential files in UPF format. The archive only contains the files and no sub directories.
 * A JSON file (`.json`) containing the recommended wavefunction and charge density cutoffs and other metadata for each pseudo potential in the archive.

## SSSP JSON file format
Each SSSP configuration comes with a file in JSON format that provides various metadata for the pseudo potential of each element in the family.
The file is a dictionary where the keys are the standard symbol for the element (one or two letters, with the first letter capitalized, e.g. "Ar" for argon) of each pseudo potential present.
For each element, the following values are provided:

    filename: the filename of the corresponding pseudo potential contained in the gzipped tarball archive.
    md5: the md5 checksum of the pseudo potential file contained in the gzipped tarball archive.
    pseudopotential: the name of the pseudo potential library from which the pseudo potential originates, e.g. "SG15"
    cutoff_wfc: the recommended cutoff for the wave functions in Rydberg
    cutoff_rho: the recommended cutoff for the charge density in Rydberg
