#!/bin/bash
# Changelog validation script based on Keep a Changelog standard

set -euo pipefail

CHANGELOG_FILE="${1:-CHANGELOG.md}"

if [ ! -f "$CHANGELOG_FILE" ]; then
    echo "❌ Changelog file not found: $CHANGELOG_FILE"
    exit 1
fi

echo "Validating changelog: $CHANGELOG_FILE"

# Check for required header elements
check_header() {
    local file="$1"
    
    if ! grep -q "# Changelog" "$file"; then
        echo "❌ Missing '# Changelog' header"
        return 1
    fi
    
    if ! grep -q "Keep a Changelog" "$file"; then
        echo "❌ Missing reference to Keep a Changelog standard"
        return 1
    fi
    
    if ! grep -q "Semantic Versioning" "$file"; then
        echo "❌ Missing reference to Semantic Versioning"
        return 1
    fi
    
    echo "✅ Header elements present"
    return 0
}

# Check for proper version format
check_versions() {
    local file="$1"
    local errors=0
    
    # Extract version headers
    grep -E "^## \[" "$file" | while read -r line; do
        version=$(echo "$line" | sed 's/## \[\(.*\)\].*/\1/')
        
        case "$version" in
            "Unreleased")
                echo "✅ Found Unreleased section"
                ;;
            *[0-9]*.[0-9]*.[0-9]*)
                if [[ $version =~ ^[0-9]+\.[0-9]+\.[0-9]+(-[a-zA-Z0-9.-]+)?$ ]]; then
                    echo "✅ Valid version format: $version"
                else
                    echo "❌ Invalid version format: $version"
                    ((errors++))
                fi
                ;;
            *)
                echo "❌ Unknown version format: $version"
                ((errors++))
                ;;
        esac
    done
    
    return $errors
}

# Check for required change types
check_change_types() {
    local file="$1"
    local sections=$(grep -E "^### (Added|Changed|Deprecated|Removed|Fixed|Security)" "$file" | wc -l)
    
    if [ "$sections" -eq 0 ]; then
        echo "❌ No standard change type sections found"
        echo "   Expected: Added, Changed, Deprecated, Removed, Fixed, Security"
        return 1
    fi
    
    echo "✅ Found $sections change type sections"
    
    # List found sections
    grep -E "^### (Added|Changed|Deprecated|Removed|Fixed|Security)" "$file" | while read -r line; do
        section=$(echo "$line" | sed 's/### //')
        echo "  - $section"
    done
    
    return 0
}

# Check for unreleased section
check_unreleased() {
    local file="$1"
    
    if ! grep -q "## \[Unreleased\]" "$file"; then
        echo "⚠️  Warning: No [Unreleased] section found"
        echo "   Consider adding one for upcoming changes"
        return 1
    fi
    
    echo "✅ Unreleased section present"
    return 0
}

# Check for proper date format
check_dates() {
    local file="$1"
    local errors=0
    
    grep -E "^## \[.*\] - [0-9]" "$file" | while read -r line; do
        date_part=$(echo "$line" | sed 's/.*- //')
        
        if [[ $date_part =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]]; then
            echo "✅ Valid date format: $date_part"
        else
            echo "❌ Invalid date format: $date_part (expected YYYY-MM-DD)"
            ((errors++))
        fi
    done
    
    return $errors
}

# Check for empty sections
check_empty_sections() {
    local file="$1"
    local in_section=false
    local section_name=""
    local section_empty=true
    local warnings=0
    
    while IFS= read -r line; do
        if [[ $line =~ ^###[[:space:]] ]]; then
            # Previous section check
            if [ "$in_section" = true ] && [ "$section_empty" = true ]; then
                echo "⚠️  Warning: Empty section: $section_name"
                ((warnings++))
            fi
            
            # New section
            section_name=$(echo "$line" | sed 's/### //')
            in_section=true
            section_empty=true
        elif [ "$in_section" = true ] && [[ $line =~ ^-[[:space:]] ]]; then
            section_empty=false
        fi
    done < "$file"
    
    # Check last section
    if [ "$in_section" = true ] && [ "$section_empty" = true ]; then
        echo "⚠️  Warning: Empty section: $section_name"
        ((warnings++))
    fi
    
    if [ $warnings -eq 0 ]; then
        echo "✅ No empty sections found"
    fi
    
    return 0
}

# Check for consistent formatting
check_formatting() {
    local file="$1"
    local errors=0
    
    # Check for proper bullet points
    if grep -E "^[[:space:]]*\*[[:space:]]" "$file" > /dev/null; then
        echo "⚠️  Warning: Found asterisk (*) bullet points, consider using dashes (-)"
    fi
    
    # Check for proper list indentation
    invalid_bullets=$(grep -E "^[[:space:]]*-[^[:space:]]" "$file" | wc -l)
    if [ "$invalid_bullets" -gt 0 ]; then
        echo "⚠️  Warning: Found $invalid_bullets bullet points without space after dash"
    fi
    
    echo "✅ Basic formatting checks passed"
    return $errors
}

# Main validation
main() {
    local file="$CHANGELOG_FILE"
    local errors=0
    
    echo "🔍 Validating changelog against Keep a Changelog standard..."
    echo ""
    
    check_header "$file" || ((errors++))
    check_versions "$file" || ((errors++))
    check_change_types "$file" || ((errors++))
    check_unreleased "$file" || true  # Warning only
    check_dates "$file" || ((errors++))
    check_empty_sections "$file" || true  # Warning only
    check_formatting "$file" || ((errors++))
    
    echo ""
    if [ $errors -eq 0 ]; then
        echo "✅ Changelog follows Keep a Changelog standard"
        return 0
    else
        echo "❌ Changelog validation failed with $errors errors"
        echo ""
        echo "For more information, see: https://keepachangelog.com/en/1.0.0/"
        return 1
    fi
}

# Show usage if requested
if [ "${1:-}" = "--help" ] || [ "${1:-}" = "-h" ]; then
    echo "Usage: $0 [CHANGELOG_FILE]"
    echo ""
    echo "Validates changelog format according to Keep a Changelog standard."
    echo "Defaults to CHANGELOG.md in current directory."
    echo ""
    echo "Checks performed:"
    echo "  - Header contains required elements"
    echo "  - Version numbers follow semantic versioning"
    echo "  - Change types follow standard categories"
    echo "  - Dates are in YYYY-MM-DD format"
    echo "  - No empty sections"
    echo "  - Consistent formatting"
    exit 0
fi

main