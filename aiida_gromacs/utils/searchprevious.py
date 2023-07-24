
from aiida import orm



def get_prev_inputs(inputs, input_labels):
    """check if input files are output files in previous processes
    and if so, add previous nodes to the inputs"""
    qb = orm.QueryBuilder()
    qb.append(orm.ProcessNode, tag='process')
    qb.order_by({orm.ProcessNode: {"ctime": "desc"}})
    added_files = []
    if qb.count() > 0:
        for entry in qb.all():
            previous_calculation = entry[0]
            for label in previous_calculation.outputs:
                if label in input_labels and label not in added_files:
                    added_files.append(label)
                    previous_output_node = previous_calculation.outputs[f"{label}"]
                    inputs[label] = previous_output_node
    return inputs