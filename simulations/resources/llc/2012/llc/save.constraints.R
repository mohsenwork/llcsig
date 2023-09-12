# Function that takes the output Z and E from enforce.rule1, enforce.rule2, or enforce.rule3
# and saves which elements of B and which elements of Ce are 0 by faithfulness, along with the 
# rule each constraint stems from. Returns a list of dictionaries, one dictionary for each constraint.
# Each dictionary should have the element i, j, rule and type.
#
# Args:
#   Z: (nxn) matrix recording which elements of B should be zero, an entry Z[i,j] > 0 means B[i,j] = 0 by faithfulness.
#   E: (nxn) matrix recording which elements of Ce should be zero, an entry E[i,j] > 0 means Ce[i,j] = 0 by faithfulness.
#   rule: Integer indicating which rule the constraints stem from (1, 2, or 3).
#   e: list specifying the experiment.
#
# Returns:
#   A list of dictionaries, one dictionary for each constraint. Each dictionary should have the element i, j, rule, type and experiment.
save.constraints <- function(Z, E, rule, e) {
  constraints <- list()
  n <- nrow(Z)
  
  # Iterate over all pairs of variables i and j
  for (i in 1:n) {
    for (j in 1:n) {

      # If Z[i,j] > 0, add a constraint to the list of constraints
      if (Z[i, j] > 0) {
        constraint <- list(i = i, j = j, rule = rule, type = 'B', e=e)
        constraints[[length(constraints) + 1]] <- constraint
      }

      # If E[i,j] > 0, add a constraint to the list of constraints
      if (E[i, j] > 0) {
        constraint <- list(i = i, j = j, rule = rule, type = 'Ce', e=e)
        constraints[[length(constraints) + 1]] <- constraint
      }  
    }
  }
  # Return the list of constraints
  return(constraints)
}
