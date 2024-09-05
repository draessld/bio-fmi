#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Get the absolute path to the directory containing this script
BASE_DIR=$(pwd)

echo "Starting project configuration..."

# 1. Install Python dependencies
echo "Installing Python dependencies..."
pip install -r requirements.txt

# 2. Compile C++ code
echo "Compiling C++ code..."
mkdir -p "$BASE_DIR/scripts/build"
cd "$BASE_DIR/scripts/build"
cmake ..
cmake --build .
cd $BASE_DIR

# 3. Run C++ tests
# echo "Compiling and running C++ tests..."
# g++ -std=c++17 -o "$BASE_DIR/scripts/build/tests" tests/cpp_tests.cpp
# "$BASE_DIR/scripts/build/tests"

# 4. Run Python tests
# echo "Running Python tests..."
# python -m unittest discover -s tests -p '*_test.py'

# 5. Create alias for run.py
echo "Creating alias 'katka2' for run.py..."
alias katka2="python $BASE_DIR/run.py"

# 6. Generate configure.json file
echo "Generating configure.json..."
cat <<EOL > configure.json
{
    "output_dir": "$BASE_DIR/output",
    "build_dir": "$BASE_DIR/scripts/build",
    "cpp_exe":{
        "build_exe": "$BASE_DIR/scripts/build/build",
        "find_exe": "$BASE_DIR/scripts/build/find",
        "kernelize_exe": "$BASE_DIR/scripts/build/kernelize",
        "minimize_exe": "$BASE_DIR/scripts/build/minimize"
    }
}
EOL

ALIAS_NAME="katka2"
ALIAS_COMMAND="python3 $BASE_DIR/run.py"

# File to modify
BASH_ALIASES_FILE="$HOME/.bash_aliases"

if [ ! -f "$BASH_ALIASES_FILE" ]; then
  touch "$BASH_ALIASES_FILE"
  echo "$BASH_ALIASES_FILE created."
fi

# Check if the alias already exists in ~/.bash_aliases
if grep -q "alias $ALIAS_NAME=" "$BASH_ALIASES_FILE" 2>/dev/null; then
  # If found, replace the existing alias
  sed -i "s|alias $ALIAS_NAME=.*|alias $ALIAS_NAME='$ALIAS_COMMAND'|" "$BASH_ALIASES_FILE"
  echo "Alias '$ALIAS_NAME' replaced in $BASH_ALIASES_FILE"
  unalias $ALIAS_NAME
  alias $ALIAS_NAME="$ALIAS_COMMAND"
else
  # If not found, add the new alias
  echo "alias $ALIAS_NAME='$ALIAS_COMMAND'" >> "$BASH_ALIASES_FILE"
  echo "Alias '$ALIAS_NAME' added to $BASH_ALIASES_FILE"
fi

source ~/.bashrc