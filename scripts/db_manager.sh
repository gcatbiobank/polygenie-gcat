#!/bin/bash

# Global paths to Python scripts
SCRIPT="pipeline/main.py"
LOG_FILE="pipeline/logs/update_db_output.log"

# Function to display the menu
show_menu() {
    printf "\n\n =====================================\n"
    printf "\tDatabase Update Manager\n"
    printf " =====================================\n\n"
    printf " 1) Update EHR Data\n"
    printf " 2) Update Individuals Data\n"
    printf " 3) Update Pipeline Data\n"
    printf " 4) Go Back\n"
    printf "\nChoose an option: "
}

# Spinner function to show a loading animation while keeping the "Processing" word
spinner() {
    local pid=$!
    local delay=0.1
    local spinstr='|/-\'
    printf "Processing "
    while [ -d /proc/$pid ]; do
        local temp=${spinstr#?}
        printf "%c" "$spinstr"
        spinstr=$temp${spinstr%"$temp"}
        sleep $delay
        printf "\b"
    done
    printf " \n"  # Ensures the spinner line ends cleanly
}

# Function to handle the menu input
handle_choice() {
    local choice
    read -r choice
    case "$choice" in
        1) run_ehr_update ;;
        2) run_individuals_update ;;
        3) run_pipeline_update ;;
        4) printf "Bye!\n"; erase_menu_content; return 1 ;;
        *) printf "Invalid option. Please try again.\n" ;;
    esac
    return 0
}

# Function to erase the menu content using a loop
erase_menu_content() {
    local i
    for (( i = 0; i < 13; i++ )); do  # Adjust the loop to clear enough lines
        printf "\033[2K"    # Clear the current line
        printf "\033[1A"    # Move cursor up one line
    done
    printf "\033[2K"        # Clear the final line after moving up
}

# Function to run the EHR update with spinner
run_ehr_update() {
    printf "Updating EHR data...\n"
    python3 "$SCRIPT" 1 0 0 0&> "$LOG_FILE" &
    spinner
    if wait $!; then
        printf "EHR data updated successfully.\n"
    else
        printf "Error: Failed to update EHR data.\n" >&2
    fi
}

# Function to run the Individuals update with spinner
run_individuals_update() {
    printf "Updating Individuals data...\n"
    python3 "$SCRIPT" 0 1 0 0&> "$LOG_FILE" &
    spinner
    if wait $!; then
        printf "Individuals data updated successfully.\n"
    else
        printf "Error: Failed to update Individuals data.\n" >&2
    fi
}

# Function to run the Pipeline update with spinner
run_pipeline_update() {
    printf "Updating Pipeline data...\n"
    python3 "$SCRIPT" 0 0 1 1&> "$LOG_FILE" &
    spinner
    if wait $!; then
        printf "Pipeline data updated successfully.\n"
    else
        printf "Error: Failed to update Pipeline data.\n" >&2
    fi
}

# Main function to drive the script
main() {
    while true; do
        show_menu
        if ! handle_choice; then
            break
        fi
    done
}

# Execute the main function
main
