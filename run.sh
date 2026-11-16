#!/bin/bash

# Colors for pretty formatting
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
MAGENTA='\033[1;35m'
NC='\033[0m' # No Color

echo -e "${BLUE}╔════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   🚀 Running Level 6 Solutions & Scores  ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════╝${NC}"

# Initialize score variables
score_2=0
score_3=0
score_4=0

# Loop through the 4 input files (1-small to 4-enormous)
for input_file in in_level-6_[1-4]*.txt; do
    # Automatically generate the matching output filename (e.g., out_level-6_1-small.txt)
    output_file="${input_file/in_/out_}"
    
    echo -e "${CYAN}▶ Processing:${NC} $input_file"
    echo -e "  ${YELLOW}➔ Output:${NC}   $output_file"
    
    # Track execution time
    start_time=$(date +%s.%N)
    
    # Run the solver
    python3 solve6b.py "$input_file" "$output_file"
    
    end_time=$(date +%s.%N)
    duration=$(awk -v a="$start_time" -v b="$end_time" 'BEGIN { printf "%.3f", b - a }')
    
    # Calculate score (number of lines in the output file)
    file_score=0
    if [ -f "$output_file" ]; then
        file_score=$(wc -l < "$output_file" | tr -d ' ')
    fi
    
    echo -e "  ${GREEN}✔ Finished in ${duration}s | Score (Lines): ${MAGENTA}$file_score${NC}\n"
    
    # Capture scores for levels 2, 3, and 4 based on filename pattern
    if [[ "$input_file" == *"_2-"* ]]; then
        score_2=$file_score
    elif [[ "$input_file" == *"_3-"* ]]; then
        score_3=$file_score
    elif [[ "$input_file" == *"_4-"* ]]; then
        score_4=$file_score
    fi
done

# Calculate the sum of levels 2, 3, and 4
total_sum=$((score_2 + score_3 + score_4))

echo -e "${BLUE}════════════════════════════════════════${NC}"
echo -e "${GREEN}✨ All output files generated successfully! ✨${NC}"
echo -e "${BLUE}════════════════════════════════════════${NC}"
echo -e "${MAGENTA}📊 Summary of Scores (Line Counts):${NC}"
echo -e "  • Level 2 (Medium):  $score_2"
echo -e "  • Level 3 (Large):    $score_3"
echo -e "  • Level 4 (Enormous): $score_4"
echo -e "${BLUE}────────────────────────────────────────${NC}"
echo -e "  ${YELLOW}🔥 Total Score (Levels 2 + 3 + 4): ${GREEN}$total_sum${NC}"
echo -e "${BLUE}════════════════════════════════════════${NC}"