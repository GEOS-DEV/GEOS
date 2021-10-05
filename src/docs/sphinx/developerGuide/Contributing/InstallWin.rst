.. _InstallWin:

Installing GEOSX on windows machine using Docker
===================================================

In this section, we will install GEOSX on a windows machine thanks to _docker_ container with precompiled version of
GEOSX's third party libraries (TPL). These steps are an adaptation of `ref:UsingDocker` for windows environement.

The different version of the containers can be found on `docker hub <https://hub.docker.com/r/geosx/ubuntu18.04-gcc8/tags?page=1&ordering=last_updated>`_.
Note that to work on the latest version of the code, you'll need to use the most recently pushed conatiner.


1. Install docker desktop
-------------------------

On your windows machine following the `steps <https://docs.docker.com/desktop/windows/install/>`_.
Download the most recent installer for Docker Desktop. Before installation let us check the current status of Windows
Subsystem Linux (WSL) on your machine as Docker will use *WSL2* as a backend. To do that, open a *PowerShell(Admin)*

.. PS admin

.. code:: shell

        PS > wsl --install
        PS > wsl --status
        PS > wsl --set-default-version 2
        PS > wsl --status

The first command should be installing *WSL2*, downloading an Uduntu distribution for it and ask for a restart.
The following commands are used to check the status, and if the *WSL* is still the default one, change it to *WSL2*.
More details on the installation procedure `here <https://docs.microsoft.com/en-us/windows/wsl/install>`_.

Once the *WSL2* is set as default, proceed with the *Docker Desktop* isntallation.

2. Start docker desktop
-------------------------

When launching Docker Desktop for the first time, you should be prompt a message informing you that it uses *WSL2*.
Using *PowerShell*, you can check that *Docker* and *WSL2* are actually running in the background

.. code:: shell

    PS > Get-Process docker
    PS > Get-Process wsl

You should be able to see one docker process and several *wsl* processes.

3. Preparing *DockerFile*
--------------------------

Let us now prepare installation picking a destination folder and editing our *Dockerfile*

.. code:: shell

    PS > cd D:/
    PS > mkdir install-geosx-docker
    PS > cd install-geosx-docker/
    PS > notepad.exe Dockerfile

Let us edit the *Dockerfile*, which is the declerative file four out container

.. code-block:: console
        
        # Define you base image for build arguments
        ARG IMG
        ARG VERSION
        ARG ORG
        FROM ${ORG}/${IMG}:${VERSION}

        # Uninstall some packages, install others.
        # I use those for clion, but VS code would have different requirements.
        # Use yum's equivalent commands for centos/red-hat images.
        # Feel free to adapt.
        RUN apt-get update
        RUN apt-get remove --purge -y texlive graphviz
        RUN apt-get install --no-install-recommends -y openssh-server gdb rsync gdbserver ninja-build

        # You will need cmake to build GEOSX.
        ARG CMAKE_VERSION=3.16.8
        RUN apt-get install -y --no-install-recommends curl ca-certificates && \
            curl -fsSL https://cmake.org/files/v${CMAKE_VERSION%.[0-9]*}/cmake-${CMAKE_VERSION}-Linux-x86_64.tar.gz | tar --directory=/usr/local --strip-components=1 -xzf - && \
            apt-get purge --auto-remove -y curl ca-certificates
        RUN apt-get autoremove -y

        # You'll most likely need ssh/sshd too (e.g. CLion and VSCode allow remote dev through ssh).
        # This is the part where I configure sshd.

        # I'm developing in a version of docker that requires root.
        # So by default I use root. Feel free to adapt.
        RUN echo "PermitRootLogin prohibit-password" >> /etc/ssh/sshd_config
        RUN echo "PermitUserEnvironment yes" >> /etc/ssh/sshd_config
        RUN mkdir -p -m 700 /root/.ssh
        # Put your own public key here!
        RUN echo "ssh-rsa [#a public ssh key]" > /root/.ssh/authorized_keys

        # Some important variables are provided through the environment.
        # You need to explicitly tell sshd to forward them.
        # Using these variables and not paths will let you adapt to different installation locations in different containers.
        # Feel free to adapt to your own convenience.
        RUN touch /root/.ssh/environment &&\
            echo "CC=${CC}" >> /root/.ssh/environment &&\
            echo "CXX=${CXX}" >> /root/.ssh/environment &&\
            echo "MPICC=${MPICC}" >> /root/.ssh/environment &&\
            echo "MPICXX=${MPICXX}" >> /root/.ssh/environment &&\
            echo "MPIEXEC=${MPIEXEC}" >> /root/.ssh/environment &&\
            echo "OMPI_CC=${CC}" >> /root/.ssh/environment &&\
            echo "OMPI_CXX=${CXX}" >> /root/.ssh/environment &&\
            echo "GEOSX_TPL_DIR=${GEOSX_TPL_DIR}" >> /root/.ssh/environment

        # This is the default ssh port that we do not need to modify.
        EXPOSE 22
        # sshd's option -D prevents it from detaching and becoming a daemon.
        # Otherwise, sshd would not block the process and `docker run` would quit.
        RUN mkdir -p /run/sshd
        ENTRYPOINT ["/usr/sbin/sshd", "-D"]

This file is pulling a docker image containing GEOSX's TPL as well as extra utils, such as *CMAKE* and preparing for ssh connexion.
Indeed, in the end, we will be able to run it in a detached mode, and connect to it to run and develop in GEOSX.

There is two things you might have noticed reading through the *Dockerfile*

- It has environment variables to be passed to it to select the proper image to pull, namely **${ORG}**, **${IMG}** and **${VERSION}**, we'll then have to declare them


    .. code:: shell

        PS> $env:VERSION='164-677'
        PS> $env:IMG='ubuntu18.04-gcc8'
        PS> $env:REMOTE_DEV_IMG="remote-dev-${env:IMG}"


Please note the preposition of *env:* in the windows formalisme. The **${ORG}** variable will be hard-coded as *geosx*. The last variable will be used as
an informative image name. Note also that the version above is the might not be the most recent version of the image at the time you are reading these lines and
should be change for the closest to the *TPL* commit associated with the GEOSX commit you cant to work from.

- You'll need to generate a ssh-key to be able to access the container without the need for defining a password. This can be done from the *PowerShell*,


    .. code:: shell

        PS > ssh-keygen.exe
        PS > cat [path-to-gen-key]/[your-key].pub


The first command will prompt message, asking for you to complete the deisred path for the key as well as a passphrase, with confirmation.
More details on `ssh-key generation <https://docs.microsoft.com/en-us/windows-server/administration/openssh/openssh_keymanagement#user-key-generation>`_.


4. Build the image and run the container
-----------------------------------------

The preliminary tasks are now done, let us build the image that will be containerized.

.. code:: shell

    PS> cd [path-to-dockerfile-folder]/
    PS > docker build --build-arg ORG=geosx --build-arg IMG=${env:IMG} --build-arg VERSION=${env:VERSION} -t ${env:REMOTE_DEV_IMG}:${env:VERSION} -f Dockerfile

As described above, we are passing our environment variables in the building stage, which offer the flwxibility of changin the version or image by a simple redefinition.
A log updating or pulling the different layer should be displayed afterwards and on the last line the *image id*. We can check that the image is created using *PowerShell* CLI:

.. code:: shell

    PS > docker images

or using the *Docker Desktop*

.. image:: win_install/win_docker_images.png
   :width: 45%

Now that we have the image build, let us run a container from,

.. code:: shell

    PS > docker run --cap-add=ALL  -d --name ${env:REMOTE_DEV_IMG}-${env:VERSION} -p 64000:22 --mount 'type=bind,source=D:/install_geosx_docker/,target=/app' ${env:REMOTE_DEV_IMG}:${env:VERSION}

Note that in addition to the detached flag (*-d*) and the name tage (*--name*), we provide *Docker* with the port the conatiner should associated to
communicate with ssh port 22, as well as a binding between a host mount point (*D:/install_geosx_docker/*) and a container mount point (*/app*) to have a peristent storage
for our development/geosx builds. More details on the `--mount options <https://docs.docker.com/storage/bind-mounts/>`_

There exact same steps can be achieve using th *Docker Desktop* GUI in the image tabs, clicking on the run button and filling the same information in the interface,

.. image:: win_install/win_docker_container.png
   :width: 45%

Coming back to our *PowerShell* terminal, we can check that our container is running and try to ssh to it.

.. code:: shell

    PS > docker ps -a
    CONTAINER ID   IMAGE                                 COMMAND               CREATED                  STATUS          PORTS                                     NAMES
    1efffac66c4c   remote-dev-ubuntu18.04-gcc8:156-642   "/usr/sbin/sshd -D"   Less than a second ago   Up 18 seconds   0.0.0.0:64000->22/tcp, :::64000->22/tcp   remote-dev-ubuntu18.04-gcc8-156-642
    PS > ssh root@localhost -p 64000
    Enter passphrase for key 'C:\************.ssh/id_rsa':
    Welcome to Ubuntu 18.04.5 LTS (GNU/Linux 5.10.16.3-microsoft-standard-WSL2 x86_64)

     * Documentation:  https://help.ubuntu.com
     * Management:     https://landscape.canonical.com
     * Support:        https://ubuntu.com/advantage
    This system has been minimized by removing packages and content that are
    not required on a system that users do not log into.

    To restore this content, you can run the 'unminimize' command.

    The programs included with the Ubuntu system are free software;
    the exact distribution terms for each program are described in the
    individual files in /usr/share/doc/*/copyright.

    Ubuntu comes with ABSOLUTELY NO WARRANTY, to the extent permitted by
    applicable law.

    root@b105f9ead860:~# cd /app && ls

We are now logged into our container and can start :ref:`QuickStart`.

.. note::
    You might be prompted that you miss certificates to clone, this can be reolved installing *ca-certificates* and updating them

        PS > apt install ca-certificates && update-ca-certificates

.. note::
    It might occur that *git-lfs* is missing then install it,

        PS > apt install git-lfs

From there you should be able to develop in your container or access it from an IDE, e.g. `VSCode <https://code.visualstudio.com/docs/remote/ssh>`_
or `MSVC19 <https://docs.microsoft.com/en-us/cpp/linux/connect-to-your-remote-linux-computer?view=msvc-160>`_.

5. Running a case
-------------------

    Once the code is configured and compiled, let us check status of the build,

.. code:: shell

    root@b105f9ead860:~# cd [path-to-build]/ && ./bin/geosx --help

Trying to launch a case using *mpirun*, you might get the following warning

.. code:: console

    root@b105f9ead860:/tmp# mpirun -np 4 /app/code/GEOSX/build-environment-debug/bin/geosx -i [geosx-case].xml -x 4 -y1 -z1
    --------------------------------------------------------------------------
    mpirun has detected an attempt to run as root.
    Running at root is *strongly* discouraged as any mistake (e.g., in
    defining TMPDIR) or bug can result in catastrophic damage to the OS
    file system, leaving your system in an unusable state.

    You can override this protection by adding the --allow-run-as-root
    option to your cmd line. However, we reiterate our strong advice
    against doing so - please do so at your own risk.
    --------------------------------------------------------------------------

A possible workaround is to create a new user and move create a run folder from this account

.. code:: shell

    root@b105f9ead860:~# adduser runner
    root@b105f9ead860:~# su runner
    runner@b105f9ead860:~# mkdir run && cd run/
    runner@b105f9ead860:~# cp [geosx-case].xml .
    runner@b105f9ead860:/tmp# mpirun -np 4 /app/code/GEOSX/build-environment-debug/bin/geosx -i [geosx-case].xml -x 4 -y 1 -z 1


.. GPU ?? https://docs.docker.com/desktop/windows/wsl/#gpu-support

