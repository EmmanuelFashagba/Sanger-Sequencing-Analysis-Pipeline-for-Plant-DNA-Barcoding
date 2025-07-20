#!/bin/bash

# Sanger Pipeline Wrapper Script
# This script activates the virtual environment and runs the pipeline

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
VENV_PATH="$HOME/sanger_env"
PIPELINE_SCRIPT="$SCRIPT_DIR/sanger_pipeline.py"

# Function to print colored output
print_color() {
    color=$1
    shift
    echo -e "${color}$@${NC}"
}

# Check if virtual environment exists
if [ ! -d "$VENV_PATH" ]; then
    print_color $RED "Error: Virtual environment not found at $VENV_PATH"
    print_color $YELLOW "Please create it with: python3 -m venv ~/sanger_env"
    exit 1
fi

# Check if pipeline script exists
if [ ! -f "$PIPELINE_SCRIPT" ]; then
    print_color $RED "Error: Pipeline script not found at $PIPELINE_SCRIPT"
    exit 1
fi

# Activate virtual environment
print_color $GREEN "Activating virtual environment..."
source "$VENV_PATH/bin/activate"

# Check if activation was successful
if [ "$VIRTUAL_ENV" != "$VENV_PATH" ]; then
    print_color $RED "Error: Failed to activate virtual environment"
    exit 1
fi

# Run the pipeline with all arguments passed to this script
print_color $GREEN "Running Sanger pipeline..."
python "$PIPELINE_SCRIPT" "$@"

# Capture exit code
exit_code=$?

# Deactivate virtual environment
deactivate

# Exit with the same code as the pipeline
exit $exit_code
