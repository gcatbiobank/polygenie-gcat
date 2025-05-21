#!/bin/bash

# Global paths to Python scripts
SCRIPT_DB="scripts/db_manager.sh"
SCRIPT_PIPELINE="scripts/launch_pipeline.sh"
SCRIPT_DEPLOYMENT="scripts/deploy_app.sh"
N_CORES=$(nproc)
GENOPRED_OPTION="target_pgs"

# Function to display the menu
show_menu() {
    printf "\n\n ===================================\n"
    printf "\tPOLYGENIE MAIN MENU\n"
    printf " ===================================\n\n"
    printf " 1) Launch Pipeline\n"
    printf " 2) Manage Database\n"
    printf " 3) Update PolyGenie Website\n"
    printf " 4) Exit\n"
    printf "\nChoose an option: "
}

# Function to handle the menu input
handle_choice() {
    local choice
    read -r choice
    case "$choice" in
        1) run_pipeline ;;
        2) erase_menu_content; ./$SCRIPT_DB ;;
        3) ./$SCRIPT_DEPLOYMENT ;;
        4) printf "Bye!\n"; return 1 ;;
        *) printf "Invalid option. Please try again.\n" ;;
    esac
    return 0
}

# Function to erase the menu content using a loop
erase_menu_content() {
    local i
    for (( i = 0; i < 12; i++ )); do  # Adjust the loop to clear enough lines
        printf "\033[2K"    # Clear the current line
        printf "\033[1A"    # Move cursor up one line
    done
    printf "\033[2K"        # Clear the final line after moving up
}


run_pipeline() {
    printf "Ensure that the GWAS sumstats are in the correct GenoPred folder and the PolyGenie config.ini file is correctly set.\n"
    printf "The pipeline is about to start, the calculations can take a while, do you want to continue [y/n]?"
    
    # Check the user's response
    read -r response 
    if [[ "$response" == "y" || "$response" == "Y" ]]; then
        ./$SCRIPT_PIPELINE $N_CORES $GENOPRED_OPTION
    else
        printf "\nOperation cancelled."
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
