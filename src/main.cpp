#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <omp.h>

typedef std::vector<double> col_t;
typedef std::vector<col_t> data_t;

struct Options
{
    std::string infile;
    std::string outfile;
    double dt = 1.0;
    double max_time = 0.0;
    std::string col_string;
    std::vector<std::pair<int,int>> correlations;
    size_t num_cols = 0;
    size_t num_cols_to_compute = 0;
};

[[noreturn]] void print_help_exit()
{
    std::cout 
        << "Usage: linear-correlator [-h] <scriptname>\n"
        << '\n'
        << "A basic linear correlator. Reads the following commands in <scriptname>.\n"
        << "Each command must be specified exactly once.\n"
        << "infile          Input data file\n"
        << "outfile         Output results file\n"
        << "dt              Time between each row of the input data file\n"
        << "max_time        Max time gap for correlations\n"
        << "columns         A string of characters, either 'x' (exclude the column),\n"
        << "                'c' (correlate the column normally), or 'f' (correlate the\n"
        << "                column after subtracting the mean) (e.g., xffffccc).\n"
        << "correlations    One pair of zero-indexed column indices for each correlation\n"
        << "                to calculate, notated as a,b for cross-correlating\n"
        << "                column a with later values of column b, e.g.,\n"
        << "                1,1 2,2 3,3 4,4 5,5 6,6 7,7 1,2 1,3 1,4 1,5 1,6 1,7 etc.\n"
        << '\n';
    exit(0);
}

Options process_script(std::istream& script)
{
    Options options;
    while (script)
    {
        std::string line;
        std::getline(script, line);
        if (line.empty())
            continue;

        std::istringstream command{ line };
        std::string keyword;
        while (command && keyword == "")
            command >> keyword;
        if (keyword == "")
        {
            std::cerr << "ERROR: Invalid line: " << line << '\n';
            exit(1);
        }

        if (keyword == "infile")
            command >> options.infile;
        else if (keyword == "outfile")
            command >> options.outfile;
        else if (keyword == "dt")
            command >> options.dt;
        else if (keyword == "max_time")
            command >> options.max_time;
        else if (keyword == "columns")
            command >> options.col_string;
        else if (keyword == "correlations")
        {
            options.correlations.clear();
            while (command)
            {
                std::string corr;
                command >> corr;
                if (corr.empty())
                    break;

                size_t pos = 0;
                auto col1 = stoi(corr, &pos);
                if (corr.size() <= pos + 1 || corr[pos] != ',')
                {
                    std::cerr << "ERROR: Invalid correlation pair: " << corr << '\n';
                    exit(1);
                }
                
                auto sub = corr.substr(pos + 1);
                auto col2 = stoi(sub);
                if (col1 < 0 || col2 < 0)
                {
                    std::cerr << "ERROR: Invalid correlation pair: " << corr << '\n';
                    exit(1);
                }

                options.correlations.emplace_back(col1, col2);
            }
            if (options.correlations.empty())
            {
                std::cerr << "ERROR: No correlations found\n";
                exit(1);
            }
        }
        else
        {
            std::cerr << "ERROR: Invalid keyword: " << keyword << '\n';
            exit(1);
        }
    }

    if (options.infile.empty())
    {
        std::cerr << "ERROR: infile not specified\n";
        exit(1);
    }

    if (options.outfile.empty())
    {
        std::cerr << "ERROR: outfile not specified\n";
        exit(1);
    }

    if (options.dt <= 0.0)
    {
        std::cerr << "ERROR: dt not specified\n";
        exit(1);
    }

    if (options.max_time == 0.0)
        options.max_time = 100.0 * options.dt;
    else if (options.max_time <= options.dt)
    {
        std::cerr << "ERROR: max_time not specified\n";
        exit(1);
    }

    if (options.col_string.empty())
    {
        std::cerr << "ERROR: columns not specified\n";
        exit(1);
    }
    options.num_cols = options.col_string.size();

    size_t num_c = 0, num_f = 0, num_x = 0;
    for (auto c : options.col_string)
    {
        if (c == 'c')
            ++num_c;
        else if (c == 'f')
            ++num_f;
        else if (c == 'x')
            ++num_x;
        else
        {
            std::cerr << "ERROR: Invalid column character: " << c << '\n';
            exit(1);
        }
    }
    if (num_x == options.col_string.size())
    {
        std::cerr << "ERROR: Column string contains only x's\n";
        exit(1);
    }

    options.num_cols_to_compute = num_c + num_f;

    return options;
}

size_t num_words(const std::string& line)
{
    std::istringstream test_iss(line);
    size_t ncols = 0;
    while (test_iss)
    {
        std::string word;
        test_iss >> word;
        if (word.empty() || word[0] == '#')
            break;
        ++ncols;
    }
    return ncols;
}

data_t input_file(std::istream& reader, const Options& options)
{
    auto num_cols = options.num_cols;
    std::vector<std::vector<std::string>> data_text(num_cols);

    auto pos = reader.tellg();
    while (reader)
    {
        pos = reader.tellg();
        std::string tmp;
        std::getline(reader, tmp);
        if (!tmp.empty() && tmp[0] != '\n' && tmp[0] != '#')
        {
            auto num = num_words(tmp);
            if (num_cols != num)
            {
                std::cerr << "ERROR: Number of characters in column string, "
                    << num_cols << ", not equal to number of columns in data file, "
                    << num << '\n';
                exit(1);
            }
            reader.seekg(pos);
            break;
        }
    }
    std::cout << "Read header\n";
    while (reader)
    {
        for (size_t i = 0; i < num_cols; ++i)
        {
            std::string tmp;
            reader >> tmp;
            if (tmp.empty())
                goto done_reading;

            if (options.col_string[i] != 'x')
                data_text[i].push_back(tmp);
        }
    }
done_reading:
    std::cout << "Read body\n";

    auto num_rows = data_text.back().size();
    data_t values(num_cols);
    for (size_t i = 0; i < num_cols; ++i)
        if (!data_text[i].empty())
            values[i].resize(num_rows);

    for (size_t i = 0; i < num_cols; ++i)
    {
        if (values[i].empty())
            continue;

#pragma omp parallel for
        for (size_t j = 0; j < num_rows; ++j)
            values[i][j] = std::stod(data_text[i][j]);
    }

    std::cout << "Read data from text\n";

    col_t avg(num_cols, 0.0);

    for (size_t i = 0; i < num_cols; ++i)
    {
        if (options.col_string[i] == 'f')
        {
            double* values_ptr = &values[i][0];
            auto& my_avg = avg[i];
#pragma omp target map(to: values_ptr) map(from: my_avg)
            {
            my_avg = 0.0;
#pragma omp teams distribute parallel for reduction(+:my_avg)
            for (size_t j = 0; j < num_rows; ++j)
                my_avg += values_ptr[j];
            }

            my_avg /= num_rows;
            std::cout << "Average value for column " << i << ": " << my_avg << '\n';
        }
    }
    std::cout << "Computed averages\n";

    for (size_t i = 0; i < num_cols; ++i)
    {
        if (avg[i])
        {
            double* values_ptr = &values[i][0];
            double& my_avg = avg[i];
#pragma omp target map(to: my_avg) map(values_ptr)
#pragma omp teams distribute parallel for
            for (size_t j = 0; j < num_rows; ++j)
                values_ptr[j] -= my_avg;
        }
    }
    std::cout << "Subtracted averages\n";

    return values;
}

data_t process_correlations(const data_t& values, const Options& options)
{
    //size_t num_cols = values.size();
    size_t num_rows = 0, idx = 0;
    while (num_rows == 0)
        num_rows = values[idx++].size();
        
    size_t num_corrs = options.correlations.size();
    size_t num_rows_per_corr = static_cast<size_t>(ceil(options.max_time / options.dt));

    data_t corr_values(num_corrs);
    for (auto& cv : corr_values)
        cv.resize(num_rows_per_corr);

    for (size_t i = 0; i < num_corrs; ++i)
    {
        std::cout << i << '\n';
        double* out_col = corr_values[i].data();
        const double* col1 = values[options.correlations[i].first].data();
        const double* col2 = values[options.correlations[i].second].data();
#pragma omp target map(to: col1, col2) map(from: out_col)
#pragma omp teams distribute
        for (size_t j = 0; j < num_rows_per_corr; ++j)
        {
            double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
            for (size_t k = 0; k < num_rows - j; ++k)
                sum += col1[k] * col2[k + j];

            out_col[j] = sum / (num_rows - j);
        }
    }

    return corr_values;
}

void write_results(
        std::ostream& writer,
        const std::vector<std::vector<double>>& corr_values,
        const double dt)
{
    for (size_t i = 0; i < corr_values.front().size(); ++i)
    {
        writer << dt * i;
        for (size_t j = 0; j < corr_values.size(); ++j)
            writer << ' ' << corr_values[j][i];
        writer << '\n';
    }
}

int main(int argc, char** argv)
{
    if (argc <= 1 || argc >= 3 || std::string(argv[1]) == "-h")
        print_help_exit();

    std::cout << "OMP number of devices: " << omp_get_num_devices() << '\n';
    std::cout << "OMP default device: " << omp_get_default_device() << '\n';

    std::string script_filename{ argv[1] };
    std::ifstream reader{ script_filename };
    if (!reader)
    {
        std::cerr << "Could not open script file: " << script_filename << '\n';
        exit(1);
    }
    auto options = process_script(reader);
    reader.close();
    std::cout << "Read script file\n";

    reader.open(options.infile);
    if (!reader)
    {
        std::cerr << "Could not open input file: " << options.infile << '\n';
        exit(1);
    }
    auto values = input_file(reader, options);
    reader.close();
    std::cout << "Read input file\n";

    auto corr_values = process_correlations(values, options);
    std::cout << "Processed correlations\n";

    std::ofstream writer{ options.outfile };
    if (!writer)
    {
        std::cerr << "Could not open output file: " << options.outfile << '\n';
        exit(1);
    }
    write_results(writer, corr_values, options.dt);
    std::cout << "Wrote results\n";
}
