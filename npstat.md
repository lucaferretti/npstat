
# Building and Using Docker for NPStat
## Introduction

NPStat is a tool designed for population genetics tests and estimators for pooled Next Generation Sequencing (NGS) data. This documentation provides detailed instructions on how to build and use the NPStat application using Docker, ensuring a consistent and reproducible environment across different systems.

## Prerequisites

Before proceeding, ensure that Docker is installed on your system. You can download and install Docker from the official Docker website.

### Installing Docker
For Linux systems, you can install Docker using the package manager of your distribution. For MacOS and Windows, you can download the Docker Desktop application from the official Docker website.
For example, to install Docker on Ubuntu/Debian or CentOS, you can use the following commands:
#### Ubuntu/Debian

```bash
sudo apt-get update
sudo apt-get install docker-ce
```
#### CentOS
```bash
sudo yum install docker
```
For more detailed installation instructions, refer to the official [Docker documentation](https://docs.docker.com/engine/install/).

#### MacOS
* Download Docker Desktop: Visit the [Docker Desktop for Mac page](https://docs.docker.com/desktop/install/mac-install/).
* Run the Installer: Open the downloaded .dmg file and drag Docker to the Applications folder.
* Start Docker Desktop: Open Docker from the Applications folder.
* Verify Installation: Open a Terminal and run:
```bash
docker --version
```
#### Windows
* Download Docker Desktop: Visit the [Docker Desktop for Windows page](https://docs.docker.com/desktop/install/windows-install/).
* Run the Installer: Open the downloaded .exe file and follow the installation instructions.



## Using NPStat container
Before using the container, you need to pull the image from Docker Hub. The image is available at the following location: [biotechvana/npstat](https://hub.docker.com/r/biotechvana/npstat).

### Pull the image from docker hub
```bash
docker pull biotechvana/npstat
```

tag the image, if you want to use a different name or to avoid using the long image name 
```bash
docker tag biotechvana/npstat npstat
```

Verify that the image has been successfully pulled by running the following command:

```bash
docker run --rm npstat
```



### Using the Docker Container
To run the NPStat application, use the following command:

```bash
docker run --rm npstat
```

This command will display the usage information for NPStat.

#### Running NPStat with Data:
To analyze specific data, you can mount a local directory containing your data files to the container and specify the input files. For example:
```bash
docker run --rm -v /path/to/data:/data npstat  -n sample_size -l window_length [options] FILE.pileup
```
* Replace `sample_size` with the haploid sample size.
* Replace `window_length` with the window length in bases.
* Replace `/path/to/data` with the path to your local directory containing the data files and `input_file` with the name of your input file. 
* **Note:** that `:/data` specifies the mount point inside the container where the data files will be available and mounted as the current working directory. And the mount point `/data` is fixed inside the container and should not be changed. `/data` is also the default working directory inside the container.

For example, if you have a data file named `FILE.pileup` in a directory `/home/user/npstat_data`, you can run NPStat as follows:
```bash
docker run --rm -v /home/user/npstat_data:/data npstat -i /data/FILE.pileup

## as /data is the working directory inside the container, you can also run it as follows
docker run --rm -v /home/user/npstat_data:/data npstat -i FILE.pileup
```
another example if you have a data file named `FILE.pileup` in your current directory, you can run NPStat as follows:
```bash
docker run --rm -v $(pwd):/data npstat -i FILE.pileup
## or
docker run --rm -v ./:/data npstat -i FILE.pileup
```

**Output:** The output file will be generated in the same directory where the input file is located.


## Building the Docker Image [Optional]
If you want to build the Docker image locally, you can follow these steps:

1. Clone the NPStat repository from GitHub:
```bash
git clone https://github.com/ahmedihafez/npstat.git
```
2. Change to the NPStat directory:
```bash
cd npstat
```
3. Build the Docker image using the provided Dockerfile:
```bash
docker build -t npstat .
```
4. Verify that the image has been successfully built by running the following command:
```bash
docker run --rm npstat
```



## Conclusion

Using Docker to build and run NPStat provides a consistent and reproducible environment, minimizing issues related to dependencies and system configurations. For more information on NPStat and its applications, refer to the NPStat GitHub repository and related literature on population genetics and NGS data analysis​ [GitHub](​https://github.com/lucaferretti/npstat).

References

* Ferretti, L., Ramos-Onsins, S.E., & Perez-Enciso, M. (2013). Population genomics from pool sequencing. Molecular Ecology, 22(17), 3845-3863. DOI: 10.1111/mec.12522.
* GitHub - lucaferretti/npstat: Population genetics tests and estimators for pooled NGS data. Retrieved from GitHub.