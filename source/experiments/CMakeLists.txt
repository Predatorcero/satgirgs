
# Check if examples are enabled
if(NOT OPTION_BUILD_EXPERIMENTS)
    return()
endif()

# Example applications
add_subdirectory(graph-testing)
