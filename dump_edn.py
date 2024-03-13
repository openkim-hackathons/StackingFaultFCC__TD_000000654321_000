# -------------------------------------------------------------------------------
# dump_string: Dumps an array of strings to a file in edn format
# -------------------------------------------------------------------------------
def dump_string_array(EdnName, string, file_id):

    Text = """
   "{}" {{
         "source-value" [ "{}" ]
   }}""".format(
        EdnName, string
    )
    file_id.write("{}\n".format(Text))


# -------------------------------------------------------------------------------
# dump_dbl_scalar: Dumps a scalar property to a file in edn format
# -------------------------------------------------------------------------------
def dump_dbl_scalar(EdnName, source_scalar, file_id, source_unit=None):

    Text = """
   "{}" {{
         "source-value" {:.15e}\n""".format(
        EdnName, source_scalar
    )

    if source_unit is not None:
        Text = Text + '         "source-unit" "{}"\n'.format(source_unit)

    Text = Text + "   }"

    file_id.write("{}\n".format(Text))


# -------------------------------------------------------------------------------
# dump_dbl_vector: Dumps a double precision vector to a file in edn format
# -------------------------------------------------------------------------------
def dump_dbl_vector(EdnName, source_vector, file_id, source_unit=None):

    str_vector = "".join(format(element, "25.15e") for element in source_vector)
    Text = """
   "{}" {{
         "source-value" [ {} ]\n""".format(
        EdnName, str_vector
    )

    if source_unit is not None:
        Text = Text + '   "source-unit" "{}"\n '.format(source_unit)

    Text = Text + "  }"

    file_id.write("{}\n".format(Text))


# -------------------------------------------------------------------------------
# dump_dbl_matrix: Dumps a double precision matrix to a file in edn format
# -------------------------------------------------------------------------------
def dump_dbl_matrix(EdnName, source_matrix, file_id, source_unit=None):

    Text = """
   "{}" {{
         "source-value" [""".format(
        EdnName
    )

    row_Idx = 0
    for row in source_matrix:

        str_row = "".join(format(element, "25.15e") for element in row)
        if row_Idx == 0:
            Text = Text + "[ {} ]\n".format(str_row)
        elif row_Idx == len(source_matrix) - 1:
            Text = Text + " " * 25 + "[ {} ]]\n".format(str_row)
        else:
            Text = Text + " " * 25 + "[ {} ]\n".format(str_row)

        row_Idx = row_Idx + 1

    if source_unit is not None:
        Text = Text + '         "source-unit" "{}"\n '.format(source_unit)

    Text = Text + "  }"

    file_id.write("{}\n".format(Text))
