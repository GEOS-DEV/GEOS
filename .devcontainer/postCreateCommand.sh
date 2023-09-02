# Activate the codespaces configuration for vscode
/usr/bin/ln -s /workspaces/GEOS/.vscode-codespaces /workspaces/GEOS/.vscode
# Activate the appropriate submodules
git submodule init
git submodule deinit integratedTests
git submodule update
