#!/bin/bash

# Phantom Science Release Script
# Automates the release process for semantic versioning

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Script configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
VERSION_FILE="$PROJECT_ROOT/VERSION"

# Functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

usage() {
    cat << EOF
Usage: $0 <version_type> [options]

Version Types:
  patch     Increment patch version (0.1.0 -> 0.1.1)
  minor     Increment minor version (0.1.0 -> 0.2.0)
  major     Increment major version (0.1.0 -> 1.0.0)
  <version> Set specific version (e.g., 1.2.3)

Options:
  --dry-run       Show what would be done without making changes
  --skip-tests    Skip running tests before release
  --skip-build    Skip building binaries
  --help         Show this help message

Examples:
  $0 patch                    # 0.1.0 -> 0.1.1
  $0 minor                    # 0.1.0 -> 0.2.0
  $0 major                    # 0.1.0 -> 1.0.0
  $0 1.0.0                    # Set to specific version
  $0 patch --dry-run          # Preview changes
  $0 minor --skip-tests       # Release without running tests
EOF
}

# Parse command line arguments
VERSION_TYPE=""
DRY_RUN=false
SKIP_TESTS=false
SKIP_BUILD=false

while [[ $# -gt 0 ]]; do
    case $1 in
        patch|minor|major)
            VERSION_TYPE="$1"
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --skip-tests)
            SKIP_TESTS=true
            shift
            ;;
        --skip-build)
            SKIP_BUILD=true
            shift
            ;;
        --help)
            usage
            exit 0
            ;;
        [0-9]*.[0-9]*.[0-9]*)
            VERSION_TYPE="$1"
            shift
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

if [[ -z "$VERSION_TYPE" ]]; then
    log_error "Version type is required"
    usage
    exit 1
fi

# Get current version
get_current_version() {
    if [[ -f "$VERSION_FILE" ]]; then
        cat "$VERSION_FILE" | sed 's/-dev$//'
    else
        echo "0.0.0"
    fi
}

# Calculate new version
calculate_new_version() {
    local current_version=$1
    local version_type=$2
    
    IFS='.' read -r major minor patch <<< "$current_version"
    
    case $version_type in
        patch)
            ((patch++))
            ;;
        minor)
            ((minor++))
            patch=0
            ;;
        major)
            ((major++))
            minor=0
            patch=0
            ;;
        [0-9]*.[0-9]*.[0-9]*)
            echo "$version_type"
            return
            ;;
        *)
            log_error "Invalid version type: $version_type"
            exit 1
            ;;
    esac
    
    echo "$major.$minor.$patch"
}

# Validate version format
validate_version() {
    local version=$1
    if ! [[ "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
        log_error "Invalid version format: $version"
        log_error "Expected format: X.Y.Z (e.g., 1.2.3)"
        exit 1
    fi
}

# Check if working directory is clean
check_git_status() {
    if ! git diff-index --quiet HEAD --; then
        log_error "Working directory is not clean. Please commit or stash changes."
        exit 1
    fi
    
    if [[ $(git rev-parse --abbrev-ref HEAD) != "main" ]]; then
        log_warning "Not on main branch. Current branch: $(git rev-parse --abbrev-ref HEAD)"
        read -p "Continue anyway? (y/N): " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            exit 1
        fi
    fi
}

# Run tests
run_tests() {
    if [[ "$SKIP_TESTS" == true ]]; then
        log_warning "Skipping tests"
        return
    fi
    
    log_info "Running tests..."
    cd "$PROJECT_ROOT"
    
    # Run Go tests
    if ! go test ./...; then
        log_error "Tests failed"
        exit 1
    fi
    
    # Run integration tests if they exist
    if [[ -f "scripts/test-integration.sh" ]]; then
        log_info "Running integration tests..."
        if ! ./scripts/test-integration.sh; then
            log_error "Integration tests failed"
            exit 1
        fi
    fi
    
    log_success "All tests passed"
}

# Build binaries
build_binaries() {
    if [[ "$SKIP_BUILD" == true ]]; then
        log_warning "Skipping build"
        return
    fi
    
    log_info "Building binaries..."
    cd "$PROJECT_ROOT"
    
    # Build for multiple platforms
    local platforms=("linux/amd64" "darwin/amd64" "darwin/arm64" "windows/amd64")
    local build_dir="build/release"
    
    mkdir -p "$build_dir"
    
    for platform in "${platforms[@]}"; do
        IFS='/' read -r os arch <<< "$platform"
        log_info "Building for $os/$arch..."
        
        local output_name="phantom-science"
        if [[ "$os" == "windows" ]]; then
            output_name="phantom-science.exe"
        fi
        
        local output_path="$build_dir/phantom-science-$os-$arch"
        mkdir -p "$output_path"
        
        if ! GOOS=$os GOARCH=$arch go build -o "$output_path/$output_name" .; then
            log_error "Build failed for $os/$arch"
            exit 1
        fi
        
        # Create tarball
        local tarball_name="phantom-science-$os-$arch.tar.gz"
        (cd "$build_dir" && tar -czf "$tarball_name" "phantom-science-$os-$arch/")
        log_success "Created $tarball_name"
    done
    
    log_success "All binaries built successfully"
}

# Update version files
update_version_files() {
    local new_version=$1
    
    log_info "Updating version files..."
    
    # Update VERSION file
    echo "$new_version" > "$VERSION_FILE"
    
    # Update CHANGELOG.md
    local changelog="$PROJECT_ROOT/CHANGELOG.md"
    local today=$(date +%Y-%m-%d)
    
    # Replace [Unreleased] with actual version
    sed -i.bak "s/## \[Unreleased\] - v[0-9.]*-dev/## [$new_version] - $today/" "$changelog"
    rm -f "$changelog.bak"
    
    # Add new unreleased section
    local next_version=$(calculate_new_version "$new_version" "minor")
    local unreleased_section="## [Unreleased] - v$next_version-dev\n\n### Added\n\n### Changed\n\n### Fixed\n\n"
    
    # Insert after the changelog header
    awk -v section="$unreleased_section" '
        /^The format is based on/ {
            print
            print ""
            print section
            next
        }
        {print}
    ' "$changelog" > "$changelog.tmp" && mv "$changelog.tmp" "$changelog"
    
    log_success "Version files updated"
}

# Create git tag and commit
create_git_tag() {
    local version=$1
    
    log_info "Creating git tag and commit..."
    
    # Add changed files
    git add "$VERSION_FILE" "$PROJECT_ROOT/CHANGELOG.md"
    
    # Create commit
    git commit -m "chore: release v$version

This release includes:
- Version bump to v$version
- Updated changelog with release notes
- Production-ready binaries for multiple platforms

🤖 Generated with [Claude Code](https://claude.ai/code)

Co-Authored-By: Claude <noreply@anthropic.com>"
    
    # Create annotated tag
    git tag -a "v$version" -m "Release v$version

Phantom Science v$version - Professional genomics simulation platform with advanced steganographic fingerprinting capabilities.

See CHANGELOG.md for detailed release notes."
    
    log_success "Git tag v$version created"
}

# Push to remote
push_release() {
    local version=$1
    
    log_info "Pushing release to remote..."
    
    # Push commits
    git push origin main
    
    # Push tags
    git push origin "v$version"
    
    log_success "Release pushed to remote"
}

# Generate release notes
generate_release_notes() {
    local version=$1
    local release_notes_file="$PROJECT_ROOT/build/release/RELEASE_NOTES_v$version.md"
    
    log_info "Generating release notes..."
    
    mkdir -p "$(dirname "$release_notes_file")"
    
    cat > "$release_notes_file" << EOF
# Phantom Science v$version Release Notes

## Overview

Phantom Science v$version delivers enterprise-grade genomics simulation with advanced steganographic fingerprinting capabilities.

## Key Features

EOF
    
    # Extract changelog for this version
    awk "/^## \[$version\]/,/^## \[/" "$PROJECT_ROOT/CHANGELOG.md" | \
        grep -v "^## \[" | \
        sed 's/^### /## /' >> "$release_notes_file"
    
    cat >> "$release_notes_file" << EOF

## Installation

### Binary Download
Download the appropriate binary for your platform from the release assets.

### Build from Source
\`\`\`bash
git clone https://github.com/scttfrdmn/phantom-science.git
cd phantom-science
make build
\`\`\`

### Docker
\`\`\`bash
docker pull phantomscience/phantom-science:v$version
\`\`\`

## Verification

All release binaries are signed and can be verified using:
\`\`\`bash
# Verify signature (when available)
gpg --verify phantom-science-linux-amd64.tar.gz.sig phantom-science-linux-amd64.tar.gz
\`\`\`

## Support

- Documentation: https://github.com/scttfrdmn/phantom-science/wiki
- Issues: https://github.com/scttfrdmn/phantom-science/issues
- Discussions: https://github.com/scttfrdmn/phantom-science/discussions

---

Generated on $(date)
EOF
    
    log_success "Release notes generated: $release_notes_file"
}

# Main release process
main() {
    cd "$PROJECT_ROOT"
    
    log_info "Starting release process for Phantom Science"
    
    # Get current and new versions
    local current_version=$(get_current_version)
    local new_version=$(calculate_new_version "$current_version" "$VERSION_TYPE")
    
    validate_version "$new_version"
    
    log_info "Current version: $current_version"
    log_info "New version: $new_version"
    
    if [[ "$DRY_RUN" == true ]]; then
        log_warning "DRY RUN MODE - No changes will be made"
        log_info "Would update version from $current_version to $new_version"
        log_info "Would run tests and build binaries"
        log_info "Would create git tag v$new_version"
        log_info "Would push to remote repository"
        exit 0
    fi
    
    # Confirm release
    echo
    log_warning "This will create release v$new_version"
    read -p "Continue? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        log_info "Release cancelled"
        exit 0
    fi
    
    # Execute release steps
    check_git_status
    run_tests
    build_binaries
    update_version_files "$new_version"
    generate_release_notes "$new_version"
    create_git_tag "$new_version"
    push_release "$new_version"
    
    echo
    log_success "🎉 Release v$new_version completed successfully!"
    log_info "Next steps:"
    log_info "1. Create GitHub release at: https://github.com/scttfrdmn/phantom-science/releases"
    log_info "2. Upload release binaries from build/release/"
    log_info "3. Publish release notes"
    log_info "4. Update documentation"
}

# Run main function
main "$@"