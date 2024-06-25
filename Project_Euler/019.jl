using Dates

function count_first_sundays(start_year, end_year)
    sundays = 0
    for year in start_year:end_year
        for month in 1:12
            date = Date(year, month, 1)
            if dayofweek(date) == 7  # 7 represents Sunday in Dates
                sundays += 1
            end
        end
    end
    return sundays
end

sundays_20th_century = count_first_sundays(1901, 2000)
println("Number of Sundays on the 1st of the month in the 20th century: ", sundays_20th_century)
