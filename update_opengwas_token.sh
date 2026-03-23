#!/bin/bash
# =============================================================================
# update_opengwas_token.sh
# 快速更新 OpenGWAS JWT Token
# =============================================================================
# Usage:
#   ./update_opengwas_token.sh <new_token>
#   ./update_opengwas_token.sh -f <token_file>
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOKEN_FILE="$HOME/.Renviron"

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

show_help() {
    echo "Usage: $0 <new_token>"
    echo "       $0 -f <token_file>"
    echo ""
    echo "Options:"
    echo "  <new_token>    New JWT token string"
    echo "  -f, --file    Read token from file"
    echo "  -h, --help    Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 eyJhbGciOiJSUzI1NiIs..."
    echo "  $0 -f /path/to/token.txt"
    echo ""
    echo "Token will be saved to: $TOKEN_FILE"
}

# Parse arguments
if [ $# -eq 0 ]; then
    show_help
    exit 1
fi

TOKEN=""
FROM_FILE=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -f|--file)
            FROM_FILE=true
            shift
            ;;
        *)
            if [ -z "$TOKEN" ]; then
                TOKEN="$1"
            fi
            shift
            ;;
    esac
done

# Handle file input
if [ "$FROM_FILE" = true ]; then
    if [ -z "$TOKEN" ]; then
        echo -e "${RED}Error: No file specified${NC}"
        exit 1
    fi
    if [ ! -f "$TOKEN" ]; then
        echo -e "${RED}Error: File not found: $TOKEN${NC}"
        exit 1
    fi
    TOKEN=$(cat "$TOKEN")
fi

# Validate token
if [ -z "$TOKEN" ]; then
    echo -e "${RED}Error: No token provided${NC}"
    exit 1
fi

# Backup existing .Renviron if it exists
if [ -f "$TOKEN_FILE" ]; then
    BACKUP_FILE="$TOKEN_FILE.backup.$(date +%Y%m%d_%H%M%S)"
    cp "$TOKEN_FILE" "$BACKUP_FILE"
    echo -e "${YELLOW}Backed up existing .Renviron to: $BACKUP_FILE${NC}"
fi

# Remove existing OPENGWAS_JWT lines and add new one
grep -v "OPENGWAS_JWT" "$TOKEN_FILE" 2>/dev/null > "$TOKEN_FILE.tmp" || true
echo "OPENGWAS_JWT=$TOKEN" >> "$TOKEN_FILE.tmp"
mv "$TOKEN_FILE.tmp" "$TOKEN_FILE"

# Also update current session environment variable
export OPENGWAS_JWT="$TOKEN"

echo -e "${GREEN}✓ Token updated successfully!${NC}"
echo "  Token length: ${#TOKEN} characters"
echo "  Saved to: $TOKEN_FILE"
echo ""
echo "You can now run MR scripts. The new token will be used automatically."
