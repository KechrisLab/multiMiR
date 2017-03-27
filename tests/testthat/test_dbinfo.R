

context("Database information functions")

test_that("multimir_dbInfo returns table source metadata", {
    dbInfo_rtn <- multimir_dbInfo()

    expect_message(multimir_dbInfo(url = "test.com"))
    expect_is(dbInfo_rtn, "data.frame")
    expect_output(str(dbInfo_rtn), "14 obs.")
    expect_output(str(dbInfo_rtn), "5 variables")
    expect_output(str(dbInfo_rtn), "map_name")
    expect_output(str(dbInfo_rtn), "source_name")
    expect_output(str(dbInfo_rtn), "source_version")
    expect_output(str(dbInfo_rtn), "source_date")
    expect_output(str(dbInfo_rtn), "source_url")
})

test_that("multimir_dbInfoVersions returns database version metadata", {
    dbInfoVers_rtn <- multimir_dbInfoVersions()

    expect_message(multimir_dbInfoVersions(url = "test.com"))
    expect_is(dbInfoVers_rtn, "data.frame")
    expect_output(str(dbInfoVers_rtn), "7 variables")
    expect_output(str(dbInfoVers_rtn), "VERSION")
    expect_output(str(dbInfoVers_rtn), "UPDATED")
    expect_output(str(dbInfoVers_rtn), "RDA")
    expect_output(str(dbInfoVers_rtn), "DBNAME")
    expect_output(str(dbInfoVers_rtn), "SCHEMA")
    expect_output(str(dbInfoVers_rtn), "PUBLIC")
    expect_output(str(dbInfoVers_rtn), "TABLES")
})

test_that("multimir_dbSchema prints sql code defining mm database to console.", {
    expect_message(multimir_dbSchema(schema.file = "test.com"))
    expect_null(multimir_dbSchema())
    expect_output(multimir_dbSchema(), "CREATE TABLE")
})


