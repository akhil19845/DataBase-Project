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
long calculate_median(std::vector<long> t)
{
    sort(t.begin(),t.end());
    int len_t = t.size();
    long med;

    if (len_t % 2 != 0)
    {
        med = t[(int)(len_t/2)];
    }
    med = (t[(int)(len_t-1)/2] + t[(int)(len_t/2)]) / 2.0;
    return med;
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
bool verify_st_bound(int points_deleted, std::vector<long> interval_list, long min_cost, long lmd_d, long interval_t)
{
    long cost = points_deleted * lmd_d;
    for(int i=points_deleted; i < interval_list.size(); i++)
    {
        cost += abs(interval_t - interval_list[i]);
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


//This function is D.P based match search to get optimal path Match, minimum cost and  length of target sequence
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
    int m = 1;

    while (m <= m_upper_bound) {
        dp[0].push_back(m * lmd_a);
        op[0].push_back(1);

        for (int i = 1; i <= n; ++i) {
            long s_m = s_0 + (m - 1) * eps_t;
            long move_res = dp[i - 1][m - 1] + abs(time[i - 1] - s_m);
            long add_res = dp[i][m - 1] + lmd_a;
            long del_res = dp[i - 1][m] + lmd_d;

            if (add_res <= move_res && add_res <= del_res) {
                dp[i].push_back(add_res);
                op[i].push_back(1);
            } else if (move_res <= add_res && move_res <= del_res) {
                dp[i].push_back(move_res);
                op[i].push_back(0);
            } else {
                dp[i].push_back(del_res);
                op[i].push_back(2);
            }
        }

        if (dp[n][m] < minimum_cost) {
            minimum_cost = dp[n][m];
            optimal_length = m;
            m_upper_bound = floor(minimum_cost / lmd_a) + n;
        }

        m += 1;
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
        interval_list.push_back(t[i] - t[i - 1]);
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
    exact_rep.min_interval_t = min_interval_t;
    exact_rep.min_starting_point = min_starting_point;
    exact_rep.optimal_length = optimal_length;
    return exact_rep;
}

// This function is used to 
long calculate_custom_interval(std::vector<long> timestamps, long interval_granularity)
{
    std::vector<long> time_diff_list;
    for (int i = 1; i < timestamps.size(); i++)
    {
        time_diff_list.push_back((long)(timestamps[i] - timestamps[i - 1]));
    }
    
    long median_time_diff = calculate_median(time_diff_list);
    long custom_interval = round(median_time_diff / interval_granularity) * interval_granularity;
    
    return custom_interval;
}

//
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
    
    for (int i = 1; i < n + 1; i++)
    {
        dp_matrix[i].push_back(i * lmd_d);
        op_matrix[i].push_back(2);
    }
    
    float optimal_length = 10e8;
    float upper_bound_m = 10e8;
    float minimum_cost = 10e8;
    int m = 1;
    
    while(m <= upper_bound_m)
    {
        dp_matrix[0].push_back(m * lmd_a);
        op_matrix[0].push_back(1);
        
        for (int i = 1; i < n+ 1; i++)
        {
            long current_point, move_result, add_result, delete_result;
            current_point = starting_point + (m - 1) * interval_t;
            move_result = dp_matrix[i - 1][m - 1] + abs(timestamps[i - 1] - current_point);
            add_result = dp_matrix[i][m - 1] + lmd_a;
            delete_result = dp_matrix[i - 1][m] + lmd_d;
            
            if (move_result <= add_result && move_result <= delete_result)
            {
                dp_matrix[i].push_back(move_result);
                op_matrix[i].push_back(0);
            }
            else if (add_result <= move_result && add_result <= delete_result)
            {
                dp_matrix[i].push_back(add_result);
                op_matrix[i].push_back(1);
            }
            else
            {
                dp_matrix[i].push_back(delete_result);
                op_matrix[i].push_back(2);
            }
        }
        
        if (dp_matrix[n][m] < minimum_cost)
        {
            minimum_cost = dp_matrix[n][m];
            optimal_length = m;
            upper_bound_m = floor(minimum_cost / lmd_a) + n;
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

//
MedianApproximation median_approximation(std::vector<long> timestamps, long lmd_a, long lmd_d, long interval_granularity) // To find the median of the given data
{
    int n = timestamps.size();
    long interval_t = calculate_custom_interval(timestamps, interval_granularity);
    int n_md = floor(n / 2);
    long s_md = calculate_median(timestamps);
    std::vector<std::vector<long>> dp_l;
    std::vector<std::vector<long>> dp_r;
    std::vector<std::vector<long>> op_l;
    std::vector<std::vector<long>> op_r;
    for (int i = 0; i < n_md + 1; i++)
    {
        dp_l.push_back(std::vector<long>());
        op_l.push_back(std::vector<long>());
        dp_r.push_back(std::vector<long>());
        op_r.push_back(std::vector<long>());
    }
    dp_l[0].push_back(0);
    op_l[0].push_back(0);
    dp_r[0].push_back(0);
    op_r[0].push_back(0);
    for (int i = 1; i < n_md + 1; i++)
    {
        dp_l[i].push_back(i * lmd_d);
        op_l[i].push_back(2);
        dp_r[i].push_back(i * lmd_d);
        op_r[i].push_back(2);
    }
    float optimal_length = 10e8;
    float m_ub = 10e8;
    float minimum_cost = 10e8;
    float m = 1;
    while (m <= m_ub)
    {
        dp_l[0].push_back(m * lmd_a);
        op_l[0].push_back(1);
        dp_r[0].push_back(m * lmd_a);
        op_r[0].push_back(1);
        for (int i = 1; i < n_md + 1; i++)
        {
            long s_m_l, s_m_r, t_i_l, t_i_r;
            if (n % 2 == 1)
            {
                s_m_l = s_md - m * interval_t;
                s_m_r = s_md + m * interval_t;
                t_i_l = timestamps[(int)((n - 1) / 2) - i];
                t_i_r = timestamps[(int)((n + 1) / 2) + (i - 1)];
            }
            else
            {
                s_m_l = s_md - (m - 0.5) * interval_t;
                s_m_r = s_md + (m - 0.5) * interval_t;
                t_i_l = timestamps[(int)(n / 2) - i];
                t_i_r = timestamps[(int)(n / 2) + i - 1];
            }
            long move_res_l = dp_l[i - 1][m - 1] + abs(t_i_l - s_m_l);
            long move_res_r = dp_r[i - 1][m - 1] + abs(t_i_r - s_m_r);
            long add_res_l = dp_l[i][m - 1] + lmd_a;
            long add_res_r = dp_r[i][m - 1] + lmd_a;
            long del_res_l = dp_l[i - 1][m] + lmd_d;
            long del_res_r = dp_r[i - 1][m] + lmd_d;
            long min_res_l = std::min({move_res_l, add_res_l, del_res_l});
            if (move_res_l == min_res_l)
            {
                dp_l[i].push_back(move_res_l);
                op_l[i].push_back(0);
            }
            else if (add_res_l == min_res_l)
            {
                dp_l[i].push_back(add_res_l);
                op_l[i].push_back(1);
            }
            else
            {
                dp_l[i].push_back(del_res_l);
                op_l[i].push_back(2);
            }
            long min_res_r = std::min({move_res_r, add_res_r, del_res_r});
            if (move_res_r == min_res_r)
            {
                dp_r[i].push_back(move_res_r);
                op_r[i].push_back(0);
            }
            else if (add_res_r == min_res_r)
            {
                dp_r[i].push_back(add_res_r);
                op_r[i].push_back(1);
            }
            else
            {
                dp_r[i].push_back(del_res_r);
                op_r[i].push_back(2);
            }
        }
        if (dp_r[n_md][m] + dp_l[n_md][m] < minimum_cost)
        {
            minimum_cost = dp_r[n_md][m] + dp_l[n_md][m];
            optimal_length = m;
            if (n % 2 == 1)
                m_ub = ((int)(floor(minimum_cost / lmd_a + n)) - 1) / 2;
            else
                m_ub = ((int)(floor(minimum_cost / lmd_a + n))) / 2;
        }
        m += 1;
    }
    long starting_point;
    // Data loss because of the type convertions
    if (n % 2 == 1)
    {
        starting_point = s_md - optimal_length * interval_t;
        m = optimal_length * 2.0 + 1.0;
    }
    else
    {
        starting_point = s_md - (optimal_length - 0.5) * interval_t;
        m = optimal_length * 2.0;
    }
    struct MedianApproximation median_approx;
    median_approx.minimum_cost = minimum_cost;
    median_approx.interval_t = interval_t;
    median_approx.starting_point = starting_point;
    median_approx.m = m;
    return median_approx;
}

MedianApproximationAll median_approximation_all(std::vector<long> timestamps, long lmd_a, long lmd_d, long interval_granularity) // Calculate the approximated dataset
{
    MedianApproximation median_approx;
    StartPointApproximation start_point_approx;
    median_approx = median_approximation(timestamps, lmd_a, lmd_d, interval_granularity);
    start_point_approx = start_point_approximation(timestamps, lmd_a, lmd_d, interval_granularity);
    MedianApproximationAll median_approx_all;
    if (median_approx.minimum_cost <= start_point_approx.minimum_cost)
    {
        median_approx_all.interval_t = median_approx.interval_t;
        median_approx_all.starting_point = median_approx.starting_point;
        median_approx_all.m = median_approx.m;
    }
    else
    {
        median_approx_all.interval_t = start_point_approx.interval_t;
        median_approx_all.starting_point = start_point_approx.starting_point;
        median_approx_all.m = start_point_approx.optimal_length;
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

float calculate_cost(std::vector<long> repair, std::vector<long> truth, long lmd_a, long lmd_d) // Calculate Cost of the given two timestamp arrays
{
    lmd_a = 5;
    lmd_d = 5;
    int n = repair.size();
    int m = truth.size();
    std::vector<std::vector<long long>> dp;
    for (int i = 0; i <= n + 1; i++)
    {
        dp.push_back(std::vector<long long>());
    }
    dp[0].push_back(0);
    lmd_a = lmd_a * (truth[1] - truth[0]);
    lmd_d = lmd_d * (truth[1] - truth[0]);
    for (int i = 1; i < n + 1; i++)
    {
        dp[i].push_back(i * lmd_d);
    }
    for (int j = 1; j < m + 1; j++)
    {
        dp[0].push_back(j * lmd_a);
        for (int i = 1; i < n + 1; i++)
        {
            long long move_cost = dp[i - 1][j - 1] + abs(repair[i - 1] - truth[j - 1]);
            long long insert_cost = dp[i][j - 1] + lmd_a;
            long long deletion_cost = dp[i - 1][j] + lmd_d;
            
            if (move_cost <= insert_cost && move_cost <= deletion_cost)
            {
                dp[i].push_back(move_cost);
            }
            else if (insert_cost <= move_cost && insert_cost <= deletion_cost)
            {
                dp[i].push_back(insert_cost);
            }
            else
            {
                dp[i].push_back(deletion_cost);
            }
        }
    }
    float overall_cost = (float)dp[n][m];
    return overall_cost;
}

float calculate_rmse(std::vector<long> truth, std::vector<long> repair) // Calculate the Root Mean Squared Error
{
    int len = std::min(truth.size(), repair.size());
    long long sum = 0;
    for (int i = 0; i < len; i++)
    {
        sum += pow(abs(truth[i] - repair[i]), 2);
    }
    float res = sqrt(sum / len);
    return res;
}

float calculate_accuracy(std::vector<long> truth, std::vector<long> fault, std::vector<long> repair) // Calcuate the auccuracy of the model
{
    int min_len;
    min_len = std::min({truth.size(), fault.size(), repair.size()});
    float error = 0;
    float cost = 0;
    float inject = 0;
    for (int i = 0; i < min_len; i++)
    {
        error += (long)abs(truth[i] - repair[i]);
    }
    for (int i = 0; i < min_len; i++)
    {
        cost += (long)abs(fault[i] - repair[i]);
    }
    for (int i = 0; i < min_len; i++)
    {
        inject += (long)abs(truth[i] - fault[i]);
    }
    if (error == 0)
    {
        return 1;
    }
    float ret = (1 - (error / (cost + inject)));
    return ret;
}

float compute_metrics(std::vector<long> repair, std::vector<long> truth, std::vector<long> fault, std::string metric_name) // Calculate the chosen metric
{
    if (metric_name == "cost")
    {
        long lmd_a = 5 * (truth[1] - truth[0]);
        long lmd_d = 5 * (truth[1] - truth[0]);
        return calculate_cost(repair, truth, lmd_a, lmd_d);
    }
    else if (metric_name == "accuracy")
    {
        return calculate_accuracy(truth, fault, repair);
    }
    else if (metric_name == "rmse")
    {
        return calculate_rmse(truth, repair);
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
        std::vector<std::string> files = {"energy", "air_quality", "pm", "syn_labdata"}; 
        int user_input;
        int filecount = 1;
        std::cout << "\nBelow are the different datasets available, select the one you would like to execute:\n" << std::endl;
        for (int i = 0; i < files.size(); ++i) {
            std::cout << i + 1 << ". " << files[i] << std::endl;
        }
        std::cout << "5. Exit from the terminal";
        std::cout << "\n\nEnter your choice Here: ";
        std::cin >> user_input;
        std::cout << std::endl;
        user_input -= 1;
        if (user_input < 0 || user_input > files.size() +1) {
            std::cout << "Invalid input" << std::endl;
            continue;
        }
        if(user_input == 5){
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
        long source_values = 0;
        long time_scale;
        Exact_Repair exact_rep;
        MedianApproximationAll median_approx_all;

        // Using the provided error data, determine the precise repaired values.
        exact_rep = exact_repair(original_seq, lmd_a, lmd_d, interval_granularity, start_point_granularity, bias_d, bias_s);

        // Determine the estimated repaired values using the provided error data.
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