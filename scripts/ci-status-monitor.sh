#!/bin/bash
# CI/CD Status Monitoring Script
# Checks GitHub Actions status and provides alerts

set -euo pipefail

GITHUB_REPO="${GITHUB_REPOSITORY:-scttfrdmn/phantom-science}"
API_BASE="https://api.github.com/repos/$GITHUB_REPO"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  --watch         Monitor continuously (every 5 minutes)"
    echo "  --json          Output in JSON format"
    echo "  --webhook URL   Send alerts to webhook URL"
    echo "  --help          Show this help message"
}

check_workflow_status() {
    local workflow_id="$1"
    local format="${2:-pretty}"
    
    echo -e "${BLUE}Checking workflow: $workflow_id${NC}"
    
    # Get latest workflow runs
    response=$(curl -s "$API_BASE/actions/workflows/$workflow_id/runs?per_page=5" 2>/dev/null || echo "{}")
    
    if [ "$format" = "json" ]; then
        echo "$response" | jq '.'
        return
    fi
    
    # Parse response
    total_count=$(echo "$response" | jq -r '.total_count // 0')
    
    if [ "$total_count" -eq 0 ]; then
        echo -e "${YELLOW}  No workflow runs found${NC}"
        return
    fi
    
    echo "$response" | jq -r '.workflow_runs[] | "\(.conclusion // "running")|\(.status)|\(.created_at)|\(.head_commit.message // "" | split("\n")[0])"' | head -5 | while IFS='|' read -r conclusion status created message; do
        case "$conclusion" in
            "success")
                status_icon="${GREEN}✅"
                ;;
            "failure")
                status_icon="${RED}❌"
                ;;
            "cancelled")
                status_icon="${YELLOW}⏹️"
                ;;
            "skipped")
                status_icon="${BLUE}⏭️"
                ;;
            "null"|"running")
                status_icon="${BLUE}🔄"
                conclusion="running"
                ;;
            *)
                status_icon="${YELLOW}❓"
                ;;
        esac
        
        formatted_date=$(date -d "$created" "+%Y-%m-%d %H:%M" 2>/dev/null || echo "$created")
        commit_msg=$(echo "$message" | cut -c1-50)
        
        echo -e "  $status_icon ${NC}$conclusion ($status) - $formatted_date - $commit_msg"
    done
}

check_overall_status() {
    local format="${1:-pretty}"
    
    if [ "$format" = "json" ]; then
        curl -s "$API_BASE/actions/runs?per_page=10" | jq '{
            total_count: .total_count,
            recent_runs: [.workflow_runs[] | {
                conclusion: .conclusion,
                status: .status,
                created_at: .created_at,
                workflow_name: .name,
                commit_message: .head_commit.message
            }]
        }'
        return
    fi
    
    echo -e "${BLUE}=== CI/CD Status Overview ===${NC}"
    echo "Repository: $GITHUB_REPO"
    echo "Timestamp: $(date)"
    echo ""
    
    # Get recent workflow runs
    response=$(curl -s "$API_BASE/actions/runs?per_page=10" 2>/dev/null || echo "{}")
    
    if [ "$(echo "$response" | jq -r '.total_count // 0')" -eq 0 ]; then
        echo -e "${YELLOW}No recent workflow runs found${NC}"
        return
    fi
    
    # Count statuses
    success_count=$(echo "$response" | jq '[.workflow_runs[] | select(.conclusion == "success")] | length')
    failure_count=$(echo "$response" | jq '[.workflow_runs[] | select(.conclusion == "failure")] | length')
    running_count=$(echo "$response" | jq '[.workflow_runs[] | select(.status == "in_progress" or .status == "queued")] | length')
    
    echo -e "${GREEN}✅ Successful: $success_count${NC}"
    echo -e "${RED}❌ Failed: $failure_count${NC}"
    echo -e "${BLUE}🔄 Running: $running_count${NC}"
    echo ""
    
    # Show recent runs
    echo "Recent runs:"
    echo "$response" | jq -r '.workflow_runs[] | "\(.conclusion // "running")|\(.name)|\(.created_at)|\(.head_commit.message // "" | split("\n")[0])"' | head -5 | while IFS='|' read -r conclusion workflow created message; do
        case "$conclusion" in
            "success") icon="${GREEN}✅" ;;
            "failure") icon="${RED}❌" ;;
            "null"|"running") icon="${BLUE}🔄"; conclusion="running" ;;
            *) icon="${YELLOW}❓" ;;
        esac
        
        formatted_date=$(date -d "$created" "+%m-%d %H:%M" 2>/dev/null || echo "$created")
        commit_msg=$(echo "$message" | cut -c1-40)
        
        echo -e "  $icon ${NC}$workflow - $formatted_date - $commit_msg"
    done
    
    # Alert on failures
    if [ "$failure_count" -gt 0 ]; then
        echo ""
        echo -e "${RED}⚠️  ALERT: $failure_count recent workflow failures detected!${NC}"
        return 1
    fi
}

send_webhook_alert() {
    local webhook_url="$1"
    local message="$2"
    
    curl -X POST "$webhook_url" \
        -H "Content-Type: application/json" \
        -d "{\"text\": \"🚨 CI/CD Alert: $message\"}" \
        2>/dev/null || echo "Failed to send webhook alert"
}

# Parse command line arguments
WATCH_MODE=false
JSON_OUTPUT=false
WEBHOOK_URL=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --watch)
            WATCH_MODE=true
            shift
            ;;
        --json)
            JSON_OUTPUT=true
            shift
            ;;
        --webhook)
            WEBHOOK_URL="$2"
            shift 2
            ;;
        --help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Set output format
if [ "$JSON_OUTPUT" = true ]; then
    OUTPUT_FORMAT="json"
else
    OUTPUT_FORMAT="pretty"
fi

# Main execution
if [ "$WATCH_MODE" = true ]; then
    echo "Starting CI/CD status monitoring (checking every 5 minutes)..."
    echo "Press Ctrl+C to stop"
    
    while true; do
        if ! check_overall_status "$OUTPUT_FORMAT"; then
            if [ -n "$WEBHOOK_URL" ]; then
                send_webhook_alert "$WEBHOOK_URL" "Workflow failures detected in $GITHUB_REPO"
            fi
        fi
        
        echo ""
        echo "Next check in 5 minutes..."
        sleep 300
    done
else
    check_overall_status "$OUTPUT_FORMAT"
fi