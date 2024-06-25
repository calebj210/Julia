function number_to_words(n)
    units = ["", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine"]
    teens = ["ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen", "sixteen", "seventeen", "eighteen", "nineteen"]
    tens = ["", "", "twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety"]

    if n == 0
        return ""
    elseif n < 10
        return units[n+1]
    elseif n < 20
        return teens[n-9]
    elseif n < 100
        return tens[div(n, 10)+1] * (n % 10 != 0 ? "-" * number_to_words(n % 10) : "")
    elseif n < 1000
        return units[div(n, 100)+1] * " hundred" * (n % 100 != 0 ? " and " * number_to_words(n % 100) : "")
    else
        return "one thousand"
    end
end

function count_letters(str)
    return length(filter(c -> isletter(c), str))
end

total_letters = sum(count_letters(number_to_words(i)) for i in 1:1000)
println("The total number of letters used is: ", total_letters)  # Output: 21124
