#!/bin/bash

# Server info
REMOTE_PATH="polygenie_app"
REMOTE_USER="udnalab"
REMOTE_HOST="192.168.3.19"
REMOTE_DIR="polygenie_app"

# Project info
PROJECT_DIR="PolyGenie"
PROD_DIR="PolyGenie/production"

# Create the production directory
cd ..
mkdir -p "$PROD_DIR"

# List of files and directories to copy
FILES_TO_COPY=(
    "app/app.py"
    "app/about.py"
    "app/__init__.py"
    "app/assets/style.css"
    "app/assets/logo.png"
    "app/assets/PolyGenie.png"
    "app/assets/favicon.ico"
    "app/assets/terms_conditions.pdf"
    "app/assets/polygenie-schema.png"
    "sqlitedb/polygenie.db"
    "sqlitedb/db_handler.py"
    "environment.yml"
)

# Copy files to production folder
for item in "${FILES_TO_COPY[@]}"; do
    # Create the target directory if it doesn't exist
    TARGET_PATH=$(dirname "$PROD_DIR/$item")
    mkdir -p "$TARGET_PATH"
        
    # Copy the file
    cp "$PROJECT_DIR/$item" "$TARGET_PATH"
done

# Copy the production folder to the server
echo "Uploading project to the server..."
scp -r "$PROD_DIR/"* "$REMOTE_USER@$REMOTE_HOST:$REMOTE_PATH/"

echo "Deployment complete"
