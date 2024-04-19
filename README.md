# IndraDemo

## Initialization Steps
```
# set up virtual environment
python3 -m venv .venv
source .venv/bin/activate
# pip install
pip install -r requirements.txt
# install pygraphviz instructions - [reference](https://pygraphviz.github.io/documentation/stable/install.html)
## MacOS
brew install pygraphviz
export GRAPHVIZ_DIR="$(brew --prefix graphviz)"
pip install pygraphviz \
    --config-settings=--global-option=build_ext \
    --config-settings=--global-option="-I$GRAPHVIZ_DIR/include" \
    --config-settings=--global-option="-L$GRAPHVIZ_DIR/lib"
# install venv into jupyter notebook
ipython kernel install --user --name=venv
See step 4 in https://www.geeksforgeeks.org/using-jupyter-notebook-in-virtual-environment/
```

