using Dates

"Format DateTime as readable string."
function format_datetime(dt::DateTime)
    return Dates.format(dt, "yyyy-mm-dd HH:MM:SS")
end

"Compute time difference between two DateTime values as d-hh:mm:ss."
function time_difference(start_time::DateTime, end_time::DateTime)
    delta_seconds = convert(Int, Dates.value(end_time - start_time)) ÷ 1000  # Convert ms → s
    days    = div(delta_seconds, 86400)
    hours   = div(delta_seconds % 86400, 3600)
    minutes = div(delta_seconds % 3600, 60)
    seconds = delta_seconds % 60
    return @sprintf("%d-%02d:%02d:%02d", days, hours, minutes, seconds)
end

"Ensure a directory exists."
function ensure_directory(path::String)
    if !isdir(path)
        mkpath(path)
    end
end