FROM ghcr.io/acts-project/ubuntu2204:master

RUN useradd -m vscode

USER vscode

ADD --chown=vscode:vscode .bashrc /home/vscode