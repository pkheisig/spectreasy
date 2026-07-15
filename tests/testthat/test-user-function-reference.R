test_that("the user function reference matches every exported signature", {
    source_reference <- testthat::test_path(
        "../../inst/extdata/spectreasy_user_function_reference.R"
    )
    installed_reference <- system.file(
        "extdata/spectreasy_user_function_reference.R",
        package = "spectreasy"
    )
    reference_file <- if (file.exists(source_reference)) {
        source_reference
    } else {
        installed_reference
    }
    testthat::skip_if_not(file.exists(reference_file))
    expressions <- parse(reference_file)
    exports <- getNamespaceExports("spectreasy")
    calls <- list()

    inspect_expression <- function(x) {
        if (is.call(x)) {
            function_name <- if (is.symbol(x[[1]])) as.character(x[[1]]) else ""
            if (function_name %in% exports) {
                arguments <- as.list(x)[-1]
                argument_names <- names(arguments)
                argument_names <- as.character(argument_names[nzchar(argument_names)])
                calls[[function_name]] <<- c(calls[[function_name]], list(argument_names))
            }
            for (i in seq_along(x)[-1]) {
                inspect_expression(x[[i]])
            }
        } else if (is.expression(x) || is.pairlist(x)) {
            for (i in seq_along(x)) {
                inspect_expression(x[[i]])
            }
        }
    }
    inspect_expression(expressions)

    expect_setequal(names(calls), exports)
    for (function_name in exports) {
        function_object <- get(function_name, envir = asNamespace("spectreasy"))
        formal_names <- as.character(setdiff(names(formals(function_object)), "..."))
        call_arguments <- calls[[function_name]]
        coverage <- vapply(
            call_arguments,
            function(x) sum(formal_names %in% x),
            integer(1)
        )
        fullest_call <- call_arguments[[which.max(coverage)]]

        expect_setequal(fullest_call, formal_names)
    }
})
