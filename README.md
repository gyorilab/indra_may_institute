# IndraDemo

## Initialization Steps
```
# set up virtual environment
python3 -m venv .venv
source .venv/bin/activate
# export rust compiler for gseapy
curl https://sh.rustup.rs -sSf | sh -s -- -y
export PATH="$PATH:$HOME/.cargo/bin"
# pip install
pip install -r requirements.txt
# install venv into jupyter notebook
ipython kernel install --user --name=venv
See step 4 in https://www.geeksforgeeks.org/using-jupyter-notebook-in-virtual-environment/
# install pygraphviz instructions
https://pygraphviz.github.io/documentation/stable/install.html
```

