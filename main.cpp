#include <iostream>
#include <stack>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <ctime>
#include <time.h>
#include <array>
#include <tuple>
#include <chrono>
#include <iomanip>

// using namespace std;
using namespace std::chrono;

struct StartPointApproximation
{
    float minimum_cost;
    long interval_t;
    long starting_point;
    float optimal_length;    
};

struct MedianApproximation
{
    float minimum_cost;
    long interval_t;
    long starting_point;
    float m;
};

struct MedianApproximationAll
{
    long interval_t;
    long starting_point;
    long m;    
};

struct MatchResult {
    std::vector<std::vector<long>> Match;
    float min_cost;
    float optimal_length;
};

struct MatchSearch {
    std::vector<std::vector<long>> Match;
    float min_cost;
    float optimal_length;
};

struct Exact_Repair
{
    long min_interval_t;
    long min_starting_point;
    float optimal_length;
};

// This function calculates the median of the time intervals
long calculate_median(const std::vector<long>& t) {
    std::vector<long> sorted_t = t; 
    std::sort(sorted_t.begin(), sorted_t.end());

    const int len_t = sorted_t.size();
    long median;

    if (len_t % 2 != 0) {
        median = sorted_t[len_t / 2];
    } else {
        median = (sorted_t[len_t / 2 - 1] + sorted_t[len_t / 2]) / 2.0;
    }

    return median;
}

//This function checks whether the value is present in the array
bool check_element_exists(const std::vector<long>& arr, long value)
{
    return std::find(arr.begin(), arr.end(), value) != arr.end();
}

//This function verifies the boundary condition
bool verify_interval_bound(long interval_t, long min_cost, std::vector<long> interval_list)
{
    long cost = 0;
    int length = interval_list.size();
    for(int i=0; i<length; i++)
    {
        cost += abs(interval_t - interval_list[i]);
    } 
    return (cost <= min_cost);
}

//This function verifies the boundary condition with deleted points
bool verify_st_bound(int points_deleted, std::vector<long> interval_list, long min_cost, long lmd_d, long interval_t) {
    long cost = points_deleted * lmd_d;
    int i = points_deleted;

    while (i < interval_list.size()) {
        cost += std::abs(interval_t - interval_list[i]);
        ++i;
    }

    return (cost < min_cost);
}
//This function is used to retrieve the path(Match) from match search
MatchResult trace_back(std::vector<std::vector<long>>& operation, const std::vector<long>& time, long s_0, long eps_t, long best_value ) {
    int n = time.size();
    std::vector<std::vector<long>> Match;
    long i = n;
    long j = best_value;

    while (i > 0 && j > 0) {
        long previous_i = i - 1;
        long previous_j = j - 1;
        if (operation[i][j] == 1) {                 // Insert
            Match.push_back({-1, previous_j});
        } else if (operation[i][j] == 0) {          // Move
             Match.push_back({previous_i, previous_j});
        } else {                                    //Delete
            Match.push_back({previous_i, -1});
        }

        i = previous_i;
        j = previous_j;
    }
    return {Match, 0.0, 0.0}; // Initialize min_cost and m_best with default values
}


//This function is D.P based match search to get optimal path Match, minimum cost and  best length of target sequence
MatchSearch match_searching(std::vector<long> time, long eps_t, long s_0, long lmd_a, long lmd_d) {
    int n = time.size();
    std::vector<std::vector<long>> dp(n + 1, std::vector<long>()); 
    std::vector<std::vector<long>> op(n + 1, std::vector<long>());

    dp[0].push_back(0);
    op[0].push_back(0);

    for (int i = 1; i <= n; ++i) {
        dp[i].push_back(i * lmd_d);
        op[i].push_back(2);
    }

    float optimal_length = 10e8;
    float m_upper_bound = 10e8;
    float minimum_cost = 10e8;
    int m_counter = 1;


    while (m_counter <= m_upper_bound) {
    const long s_m = s_0 + (m_counter - 1) * eps_t;

    dp[0].emplace_back(m_counter * lmd_a);
    op[0].emplace_back(1);

    for (int i = 1; i <= n; ++i) {
        const long diff_time_s_m = std::abs(time[i - 1] - s_m);
        const long add_res = dp[i][m_counter - 1] + lmd_a;
        const long move_res = dp[i - 1][m_counter - 1] + std::abs(time[i - 1] - s_m);
        const long del_res = dp[i - 1][m_counter] + lmd_d;

        const long min_cost = std::min({add_res, move_res, del_res});
        dp[i].emplace_back(min_cost);

        if (min_cost == add_res) {
        op[i].emplace_back(1);
        } else if (min_cost == move_res) {
         op[i].emplace_back(0);
        } else {
        op[i].emplace_back(2);
        }
    }
    const long current_cost = dp[n][m_counter];
    if (current_cost < minimum_cost) {
        optimal_length = m_counter;
        minimum_cost = current_cost;
        m_upper_bound = floor(minimum_cost / lmd_a) + n;
    }

    ++m_counter;
}
    auto Match = trace_back(op, time, s_0, eps_t, optimal_length);
    return {Match.Match, minimum_cost, optimal_length};
}

// This function is used to approximate the median value
long round_to_granularity(long interval_med, long interval_granularity)
{
    return round(interval_med / interval_granularity) * interval_granularity;
}

// This function determines optimal starting point and time interval of the target sequence
Exact_Repair exact_repair(std::vector<long> t, long lmd_a, long lmd_d, int interval_granularity, int start_point_granularity, int bias_d, int bias_s) // Exact Repair Algorithm
{
    std::vector<long> interval_list;
    int n = t.size();

    for (int i = 1; i < n - 1; i++)
    {
        int timestamp_diff = t[i] - t[i-1];
        interval_list.push_back(timestamp_diff);
        i += 1;
    }
    long interval_md = calculate_median(interval_list);
    long interval_t = round_to_granularity(interval_md, interval_granularity);
    std::vector<long> interval_t_traverse_range;
    std::vector<long> interval_t_traversed;
    float min_cost = 10e8;
    long min_interval_t;
    long min_starting_point;
    bool flag_increase = false;
    bool flag_decrease = false;
    bool min_cost_change = false;
    float optimal_length;
    int deleted_points;
    long starting_point;
    MatchSearch match_search;
    long start_ub;
    long start_lb;
    while (true)
    {
        deleted_points = 0;
        while (deleted_points < interval_list.size() && deleted_points < bias_d && ((deleted_points == 0) || verify_st_bound(deleted_points, interval_list, min_cost, lmd_d, interval_t)) )
        {
            starting_point = t[deleted_points];
            flag_increase = false;
            flag_decrease = false;
            min_cost_change = false;
            start_ub = t[deleted_points] + bias_s;
            while (starting_point <= start_ub)
            {
                match_search = match_searching(t, interval_t, starting_point, lmd_a, lmd_d);
                if (match_search.min_cost >= min_cost)
                {
                    flag_increase = true;
                    break;
                }
                else
                {
                    min_cost = match_search.min_cost;
                    optimal_length = match_search.optimal_length;
                    min_interval_t = interval_t;
                    min_starting_point = starting_point;
                    min_cost_change = true;
                    starting_point += start_point_granularity;
                }
            }
            starting_point = t[deleted_points] - 1;
            start_lb = t[deleted_points] - bias_s;
            while (starting_point >= start_lb)
            {
                starting_point -= start_point_granularity;
                match_search = match_searching(t, interval_t, starting_point, lmd_a, lmd_d);
                if(match_search.min_cost >= min_cost)
                {
                    flag_decrease = true;
                    break;                    
                }
                else
                {
                    min_cost = match_search.min_cost;
                    optimal_length = match_search.optimal_length;
                    min_interval_t = interval_t;
                    min_starting_point = starting_point;
                    min_cost_change = true;
                }
            }
            if (flag_increase && flag_decrease)
            {
                break;
            }
            deleted_points += 1;
        }
        if (verify_interval_bound(interval_t, min_cost, interval_list) == false || min_cost_change == false)
        {
            break;
        }
        
        if (check_element_exists(interval_t_traversed, (interval_t + interval_granularity)) == false && (interval_t + interval_granularity) <= round_to_granularity(interval_md, interval_granularity) + interval_granularity)
        {
            interval_t_traverse_range.push_back((interval_t + interval_granularity));
        }
        
        if (check_element_exists(interval_t_traversed, (interval_t - interval_granularity)) == false && (interval_t - interval_granularity) >= round_to_granularity(interval_md, interval_granularity) + interval_granularity)
        {
            interval_t_traverse_range.push_back((interval_t - interval_granularity));
        }

        interval_t_traversed.push_back(interval_t);
        
        if (interval_t_traverse_range.size() == 0)
        {
            break;
        }
        interval_t = interval_t_traverse_range.back();
        interval_t_traverse_range.pop_back();
    }
    Exact_Repair exact_rep;
    exact_rep.min_starting_point = min_starting_point;
    exact_rep.optimal_length = optimal_length;
    exact_rep.min_interval_t = min_interval_t;
    return exact_rep;
}


long calculate_custom_interval(std::vector<long> timestamps, long interval_granularity)
{
    std::vector<long> time_diff_list;
    int i=1;
    while(i < timestamps.size()){
    //for (int i = 1; i < timestamps.size(); i++)
   // {
        time_diff_list.push_back((long)(timestamps[i] - timestamps[i - 1]));
        i++;
    }
    
    long median_time_diff = calculate_median(time_diff_list);
    long custom_interval = round(median_time_diff / interval_granularity) * interval_granularity;
    
    return custom_interval;
}

StartPointApproximation start_point_approximation(std::vector<long> timestamps, long lmd_a, long lmd_d, long interval_granularity)
{
    long starting_point = timestamps[0];
    int n = timestamps.size();
    long interval_t = calculate_custom_interval(timestamps, interval_granularity);
    std::vector<std::vector<long>> dp_matrix;
    std::vector<std::vector<long>> op_matrix;
    
   

    for (int i = 0; i < n + 1; i++)
    {
        dp_matrix.push_back(std::vector<long>());
        op_matrix.push_back(std::vector<long>());
    }
    dp_matrix[0].push_back(0);
    op_matrix[0].push_back(0);
    

    int i=1;
    while(i < n+1){
        dp_matrix[i].push_back(i * lmd_d);
        op_matrix[i].push_back(2);
        i++;
    }
    float optimal_length = 10e8;
    float upper_bound_m = 10e8;
    float minimum_cost = 10e8;
    int m = 1;
    
    while (m <= upper_bound_m) {
    dp_matrix[0].push_back(m * lmd_a);
    op_matrix[0].push_back(1);

    int i=1;
    while(i < n + 1){
        long current_point = starting_point + (m - 1) * interval_t;
        long move_result = dp_matrix[i - 1][m - 1] + std::abs(timestamps[i - 1] - current_point);
        long add_result = dp_matrix[i][m - 1] + lmd_a;
        long delete_result = dp_matrix[i - 1][m] + lmd_d;

        int opCode;
        long minResult = std::min({move_result, add_result, delete_result});

        switch (minResult) {
            case 0:
                opCode = 0;
                break;
            case 1:
                opCode = 1;
                break;
            case 2:
                opCode = 2;
                break;
        }

        dp_matrix[i].push_back((op_matrix[i].emplace_back(opCode), minResult));
        i++;
    }

    long dpNnM = dp_matrix[n][m];
    if (dpNnM < minimum_cost) {
        minimum_cost = dpNnM;
        optimal_length = m;
        upper_bound_m = std::floor(minimum_cost / lmd_a) + n;
    }

    m += 1;
}
    
    
    struct StartPointApproximation start_point_approx;
    start_point_approx.minimum_cost = minimum_cost;
    start_point_approx.interval_t = interval_t;
    start_point_approx.starting_point = starting_point;
    start_point_approx.optimal_length= optimal_length;
    
    return start_point_approx;
}

MedianApproximation median_approximation(std::vector<long> timestamps, long lmd_a, long lmd_d, long interval_granularity) // To find the median of the given data
{
    int n = timestamps.size();
    long interval_t = calculate_custom_interval(timestamps, interval_granularity);
    int n_md = floor(n / 2);
    long s_md = calculate_median(timestamps);
    std::vector<std::vector<long>> leftdyp;
    std::vector<std::vector<long>> rightdyp;
    std::vector<std::vector<long>> leftopr;
    std::vector<std::vector<long>> rightopr;
    int i=1;
    while(i < n_md + 1){
        leftdyp.push_back(std::vector<long>());
        leftopr.push_back(std::vector<long>());
        rightdyp.push_back(std::vector<long>());
        rightopr.push_back(std::vector<long>());
        i++;
    }
    leftdyp[0].push_back(0);
    leftopr[0].push_back(0);
    rightdyp[0].push_back(0);
    rightopr[0].push_back(0);
   
    i=1;
    while(i < n_md +1){
        leftdyp[i].push_back(i * lmd_d);
        leftopr[i].push_back(2);
        rightdyp[i].push_back(i * lmd_d);
        rightopr[i].push_back(2);
        i++;
    }
    float optimal_length = 10e8;
    float m_ub = 10e8;
    float minimum_cost = 10e8;
    float m = 1;

    while (m <= m_ub) {
    long mA = m * lmd_a;
    leftdyp[0].push_back(mA);
    leftopr[0].push_back(1);
    rightdyp[0].push_back(mA);
    rightopr[0].push_back(1);

    int i=1;
    while (i < n_md + 1){
    //for (int i = 1; i < n_md + 1; i++) {
        long s_m_l, s_m_r, t_i_l, t_i_r;
        if (n % 2 != 1) {
            s_m_l = s_md - (m - 0.5) * interval_t;
            s_m_r = s_md + (m - 0.5) * interval_t;
            t_i_l = timestamps[n / 2 - i];
            t_i_r = timestamps[n / 2 + i - 1];
        } else {
            s_m_l = s_md - m * interval_t;
            s_m_r = s_md + m * interval_t;
            t_i_l = timestamps[(n - 1) / 2 - i];
            t_i_r = timestamps[(n + 1) / 2 + (i - 1)];
        }

        long move_res_l = leftdyp[i - 1][m - 1] + std::abs(t_i_l - s_m_l);
        long move_res_r = rightdyp[i - 1][m - 1] + std::abs(t_i_r - s_m_r);
        long add_res_l = leftdyp[i][m - 1] + lmd_a;
        long add_res_r = rightdyp[i][m - 1] + lmd_a;
        long del_res_l = leftdyp[i - 1][m] + lmd_d;
        long del_res_r = rightdyp[i - 1][m] + lmd_d;

        long min_res_l = std::min({move_res_l, add_res_l, del_res_l});
        long min_res_r = std::min({move_res_r, add_res_r, del_res_r});

        leftdyp[i].push_back((leftopr[i].emplace_back(move_res_l == min_res_l ? 0 : (add_res_l == min_res_l ? 1 : 2)), min_res_l));
        rightdyp[i].push_back((rightopr[i].emplace_back(move_res_r == min_res_r ? 0 : (add_res_r == min_res_r ? 1 : 2)), min_res_r));
        i++;
    }

    long dpRnMdM = rightdyp[n_md][m];
    long dpLnMdM = leftdyp[n_md][m];
    if (dpRnMdM + dpLnMdM < minimum_cost) {
        minimum_cost = dpRnMdM + dpLnMdM;
        optimal_length = m;
        int floorMinCostLaN = (int)(floor(minimum_cost / lmd_a + n));
        m_ub = (n % 2 == 1) ? (floorMinCostLaN - 1) / 2 : floorMinCostLaN / 2;
    }

    m += 1;
}
long starting_point;

constexpr double MULTIPLIER_ODD = 2.0;
constexpr double OFFSET_ODD = 1.0;
constexpr double OFFSET_EVEN = 0.5;

if (n % 2 == 1) {
    starting_point = s_md - optimal_length * interval_t;
    m = optimal_length * MULTIPLIER_ODD + OFFSET_ODD;
} else {
    starting_point = s_md - (optimal_length - OFFSET_EVEN) * interval_t;
    m = optimal_length * MULTIPLIER_ODD;
}

    struct MedianApproximation median_approx;
    median_approx.minimum_cost = minimum_cost;
    median_approx.interval_t = interval_t;
    median_approx.starting_point = starting_point;
    median_approx.m = m;
    return median_approx;
}
MedianApproximationAll median_approximation_all(std::vector<long> timestamps, long lmd_a, long lmd_d, long interval_granularity) {
    MedianApproximation median_approx = median_approximation(timestamps, lmd_a, lmd_d, interval_granularity);
    StartPointApproximation start_point_approx = start_point_approximation(timestamps, lmd_a, lmd_d, interval_granularity);

    MedianApproximationAll median_approx_all;

    if (median_approx.minimum_cost > start_point_approx.minimum_cost) {
        std::tie(median_approx_all.interval_t, median_approx_all.starting_point, median_approx_all.m) =
            std::tie(start_point_approx.interval_t, start_point_approx.starting_point, start_point_approx.optimal_length);
    } else {
        std::tie(median_approx_all.interval_t, median_approx_all.starting_point, median_approx_all.m) =
            std::tie(median_approx.interval_t, median_approx.starting_point, median_approx.m);
    }

    return median_approx_all;
}

std::vector<long> equal_series_generate(long interval_t, long starting_point, long m) // Generate Equal Series Data
{
    std::vector<long> ret;
    for (int i = 0; i < m; i++)
    {
        long value = (long)(starting_point + i * interval_t);
        ret.push_back(value);
    }
    return ret;
}

float calculate_cost(const std::vector<long>& repair, const std::vector<long>& truth, long lmd_a, long lmd_d) {
    int n = repair.size();
    int m = truth.size();
    std::vector<std::vector<long long>> dp(n + 1, std::vector<long long>(m + 1, 0));

    for (int i = 1; i <= n; i++) {
        dp[i][0] = i * lmd_d;
    }

    for (int j = 1; j <= m; j++) {
        dp[0][j] = j * lmd_a;
        for (int i = 1; i <= n; i++) {
            long long move_cost = dp[i - 1][j - 1] + std::abs(repair[i - 1] - truth[j - 1]);
            long long insert_cost = dp[i][j - 1] + lmd_a;
            long long deletion_cost = dp[i - 1][j] + lmd_d;

            dp[i][j] = std::min({move_cost, insert_cost, deletion_cost});
        }
    }

    return static_cast<float>(dp[n][m]);
}

float calculate_rmse(std::vector<long> truth, std::vector<long> repair) // Calculate the Root Mean Squared Error
{
    int length = std::min(truth.size(), repair.size());
    long long sum = 0;
    int i = 0;

    while (i < length) {
        sum += std::pow(std::abs(truth[i] - repair[i]), 2);
        ++i;
    }

    float res = std::sqrt(static_cast<float>(sum) / length);
    return res;
    
}

float calculate_accuracy(const std::vector<long>& truth, const std::vector<long>& fault, const std::vector<long>& repair) {
    const int min_len = std::min({truth.size(), fault.size(), repair.size()});

    float error = 0.0f, cost = 0.0f, inject = 0.0f;

    for (int i = 0; i < min_len; ++i) {
        const long diff_truth_repair = std::abs(truth[i] - repair[i]);
        const long diff_fault_repair = std::abs(fault[i] - repair[i]);
        const long diff_truth_fault = std::abs(truth[i] - fault[i]);

        error += static_cast<float>(diff_truth_repair);
        cost += static_cast<float>(diff_fault_repair);
        inject += static_cast<float>(diff_truth_fault);
    }

    return (error == 0.0f) ? 1.0f : (1.0f - (error / (cost + inject)));
}
float compute_metrics(const std::vector<long>& repair, const std::vector<long>& truth, const std::vector<long>& fault, const std::string& metric_name) {
    long lmd_a = 5 * (truth[1] - truth[0]);
    long lmd_d = 5 * (truth[1] - truth[0]);

    if (metric_name == "rmse") {
        return calculate_rmse(truth, repair);
    } else if (metric_name == "cost") {
        return calculate_cost(repair, truth, lmd_a, lmd_d);
    } else if (metric_name == "accuracy") {
        return calculate_accuracy(truth, fault, repair);
    }

    return 0;
}



// Function to print metrics in table format
void printMetric(std::string metricName, double resultExact, double resultApprox) {
    std::cout << std::left << std::setw(30) << metricName
         << std::setw(10) << "Exact Algorithm: " << resultExact
         << std::setw(20) << "\tApproximation Algorithm: " << resultApprox
         << std::endl;
}

int main() // Main Function
{
    while (true)
    {
        // Datasets
        std::vector<std::string> files = {"energy", "air_quality", "pm"};  
        int user_input;
        int filecount = 1;
        std::cout << "\nBelow are the different datasets available, select the one you would like to execute:\n" << std::endl;
        for (int i = 0; i < files.size(); ++i) {
            std::cout << i + 1 << ". " << files[i] << std::endl;
        }
        std::cout << "4. Exit from the terminal";
        std::cout << "\n\nEnter your choice Here: ";
        std::cin >> user_input;
        std::cout << std::endl;
        user_input -= 1;
        if (user_input < 0 || user_input > files.size() +1) {
            std::cout << "Invalid input" << std::endl;
            continue;
        }
        if(user_input == 4){
            exit(0);
        }
        
        std::cout << "You have chosen " << files[user_input] << "' dataset" << std::endl;
        std::cout << "Out of five files, which one would you like to execute?: ";
        std::cin >> filecount;
        std::cout << std::endl;
        std::cout << std::endl;
        filecount -= 1;
        std::cout <<"Executing data set file: "<< filecount + 1;
        std::cout << std::endl;
        auto start = high_resolution_clock::now();
        std::string line, word;
        int start_point_granularity = 1;
        int interval_granularity = 1;
        int lmd_a = 100;
        int lmd_d = 100;
        int bias_d = 1;
        int bias_s = 3;
        double result_exact_rmse = 0.0;
        double result_approx_rmse = 0.0;
        double result_exact_acc = 0.0;
        double result_approx_acc = 0.0;
        double result_exact_cost = 0.0;
        double result_approx_cost = 0.0;
        std::string file_name = "./data/" + files[user_input] + "/series_" + std::to_string(filecount) + ".csv";
        std::vector<long> original_seq;
        std::vector<long> ground_truth_seq;
        std::string metric;
        std::ifstream inputFile;
        std::ifstream inputFile1;
        std::ifstream inputFile2;
        // There are two files for Original and Ground Truth Values if Energy is selected as the input.
        if (user_input == 0)
        {
            file_name = "./data/dirty_energy/series_" + std::to_string(filecount) + ".csv";
            std::string data_truth = "./data/energy/series_" + std::to_string(filecount) + ".csv";
            inputFile1.open(file_name);
            getline(inputFile1, line);
            line = "";
            while (getline(inputFile1, line))
            {
                std::string first_column;
                std::string second_column;
                std::stringstream inputString(line);
                getline(inputString, first_column, ',');
                getline(inputString, second_column);
                original_seq.push_back(stol(second_column));
                line = "";
            }
            inputFile2.open(data_truth);
            getline(inputFile2, line);
            line = "";
            while (getline(inputFile2, line))
            {
                std::string first_column;
                std::string second_column;
                std::string thrid_column;
                std::stringstream inputString(line);
                getline(inputString, first_column, ',');
                getline(inputString, second_column, ',');
                getline(inputString, thrid_column);
                ground_truth_seq.push_back(stol(second_column));
                line = "";
            }
        }
        else // They have both the original and ground truth in a single file for every other dataset.
        {
            inputFile.open(file_name);
            getline(inputFile, line);
            line = "";
            while (getline(inputFile, line))
            {
                std::string first_column;
                std::string second_column;
                std::string third_column;
                std::stringstream inputString(line);
                getline(inputString, first_column, ',');
                getline(inputString, second_column, ',');
                getline(inputString, third_column);
                ground_truth_seq.push_back(stol(second_column));
                original_seq.push_back(stol(third_column));
                line = "";
            }
        }
        Exact_Repair exact_rep;
        MedianApproximationAll median_approx_all;

        // Find the exact repaired values using the error data that has been provided.
        exact_rep = exact_repair(original_seq, lmd_a, lmd_d, interval_granularity, start_point_granularity, bias_d, bias_s);

        // Using the supplied error data, calculate the estimated repaired values.
        median_approx_all = median_approximation_all(original_seq, lmd_a, lmd_d, interval_granularity);
        
        //In order to produce the final data for evaluation
        std::cout << "\n-------------------Generating Values-------------------\n" << std::endl;
        std::vector<long> exact_res = equal_series_generate(exact_rep.min_interval_t, exact_rep.min_starting_point, exact_rep.optimal_length);
        std::vector<long> appro_res = equal_series_generate(median_approx_all.interval_t, median_approx_all.starting_point, median_approx_all.m);
        

        // Calculate metrics
         result_exact_rmse = compute_metrics(exact_res, ground_truth_seq, original_seq, "rmse");
         result_approx_rmse = compute_metrics(appro_res, ground_truth_seq, original_seq, "rmse");

         result_exact_acc = compute_metrics(exact_res, ground_truth_seq, original_seq, "accuracy");
         result_approx_acc = compute_metrics(appro_res, ground_truth_seq, original_seq, "accuracy");

         result_exact_cost = compute_metrics(exact_res, ground_truth_seq, original_seq, "cost");
         result_approx_cost = compute_metrics(appro_res, ground_truth_seq, original_seq, "cost");

    // Print metrics in table format
    printMetric("a. RMSE Metrics:", result_exact_rmse, result_approx_rmse);
    printMetric("b. Accuracy:", result_exact_acc, result_approx_acc);
    printMetric("c. Total Cost:", result_exact_cost, result_approx_cost);

    // Calculate execution time
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    // Print execution time
    std::cout << std::left << std::setw(20) << "\nTotal Time taken for execution:" << duration.count() << " microseconds\n" << std::endl;
        
    }
    return 0;
}