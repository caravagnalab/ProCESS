/*
 * This file is part of the ProCESS (https://github.com/caravagnalab/ProCESS/).
 * Copyright (c) 2023-2026 Alberto Casagrande <alberto.casagrande@uniud.it>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <vector>

#include <Rcpp.h>

#include "logics_impl.hpp"

using namespace Rcpp;

RCPP_MODULE(Logics)
{

//' @name Variable
//' @title Represent a simulation quantity
//' @description The objects of this class represent one among the following
//'   quantities:
//'
//'   -  the cardinality of a species;
//'   -  the number fired event among deaths, duplications and switches in a
//'   species;
//'   -  the elapse simulation time.
//'
//' @examples
//' # build a simulation and add two species to it
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation with epigenetic states
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add the mutant "A" and set its species rates
//' sim$add_mutant("A", list(E1 = list(duplication = 0.2, death = 0.1,
//'                                    E2 = 0.01),
//'                          E2 = list(duplication = 0.08, death = 0.01,
//'                                    E1 = 0.01)))
//'
//' # get the variable representing the cardinality of "A[E1]" in sim
//' sim$var("A[E1]")
//'
//' # get the variable representing the number of duplications in "A[E1]"
//' sim$var("A[E1].duplications")
//'
//' # get the variable representing the simulation time
//' sim$var("Time")
//'
//' # the logic variables can be stored in an R variable
//' va_p <- sim$var("A[E1]")
//' va_p
//' @seealso [TissueSimulation$var()] to build a variable
    class_<Logics::Variable>("Variable")
        .method("show", &Logics::Variable::show, "Show a variable");

//' @name Expression
//' @title Represent an expression of value and variable
//' @description The objects of this class are polynomial expression of
//'   variables and values built by using the binary operators `+`, `-`,
//'   and `*` with the usual semantics.
//'
//'   An expression is one of the following object:
//'
//'   -  a variable;
//'   -  a numeric value, e.g., `3.4`;
//'   -  the sum of two expressions;
//'   -  the subtraction of two expressions;
//'   -  the multiplication of two expressions.
//'
//' @examples
//' # build a simulation and add two species to it
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation
//' sim <- TissueSimulation()
//' sim$add_mutant("A", c(duplication = 0.2, death = 0.1))
//'
//' # build an expression
//' sim$var("A") - 2 * sim$var("Time") * sim$var("A.duplications") + 3.4
//'
//' # R variables storing logic variables can also be used in expressions
//' v_time <- sim$var("Time")
//' sim$var("A") - 2 * v_time * sim$var("A.duplications") + 3.4
//'
//' # the logic expression can be stored in an R variable
//' v_exp <- sim$var("A") - 2 * v_time * sim$var("A.duplications") + 3.4
//' v_exp
//' @seealso [Variable()]
    class_<Logics::Expression>("Expression")
        .method("show", &Logics::Expression::show, "Show an expression");

    function("logics_sum", (SEXP (*)(const SEXP &, const SEXP &))&Logics::sum,
             "Sums two expressions");
    function("logics_subtract", (SEXP (*)(const SEXP &, const SEXP &))&Logics::subtract,
             "Sums two expressions");
    function("logics_multiply", (SEXP (*)(const SEXP &, const SEXP &))&Logics::multiply,
             "Multiplies two expressions");

    function("logics_lt", (SEXP (*)(const SEXP &, const SEXP &))&Logics::lt,
             "Builds a \"less than\" relation");
    function("logics_le", (SEXP (*)(const SEXP &, const SEXP &))&Logics::le,
             "Builds a \"less than or equal to\" relation");
    function("logics_eq", (SEXP (*)(const SEXP &, const SEXP &))&Logics::eq,
             "Builds an \"equal to\" relation");
    function("logics_ne", (SEXP (*)(const SEXP &, const SEXP &))&Logics::ne,
             "Builds an \"not equal to\" relation");
    function("logics_ge", (SEXP (*)(const SEXP &, const SEXP &))&Logics::ge,
             "Builds an \"greater than or equal to\" relation");
    function("logics_gt", (SEXP (*)(const SEXP &, const SEXP &))&Logics::gt,
             "Builds an \"greater than\" relation");

//' @name Formula
//' @title First-order simulation state formulas
//' @description The objects of this class describes properties of the
//'   simulation status.
//'
//'   A formula is:
//'
//'   -  a relation among two expressions (operators `<`, `<=`, `==`,
//'   `!=`, `>=`, `>`);
//'   -  the conjunction of two formulas (operator `&`);
//'   -  the disjunction of two formulas (operator `|`).
//'
//' @examples
//' # build a simulation and add two species to it
//' # set the seed of the random number generator
//' set.seed(0)
//'
//' # create a simulation with epigenetic states
//' sim <- TissueSimulation(epigenetic_states = c("E1", "E2"))
//'
//' # add the mutant "A" and set its species rates
//' sim$add_mutant("A", list(E1 = list(duplication = 0.2, death = 0.1,
//'                                    E2 = 0.01),
//'                          E2 = list(duplication = 0.08, death = 0.01,
//'                                    E1 = 0.01)))
//'
//' # get a formula that holds when the cardinality of the mutant A
//' # is greater than 1000
//' f1 <- sim$var("A[E1]") + sim$var("A[E2]") > 1000
//'
//' # get a formula that holds when the simulated time is 10 at least
//' f2 <- sim$var("Time") >= 40
//'
//' # get a formula that holds when the number of duplications doubles
//' # the switch from A[E1]
//' f3 <- sim$var("A[E1].duplications") > 2 * sim$var("A[E1].switches")
//'
//' # combine above formulas by using Boolean operators `&` and `|`
//' f1 & (f2 | f3)
//' @seealso [Variable()], [TissueSimulation$var()], `vignette("tissue_simulation")`
    class_<Logics::Formula>("Formula").method("show", &Logics::Formula::show,
                                              "Show a formula");

    function("not", (SEXP (*)(const SEXP &))&Logics::logics_not,
             "Builds the negation of a formula");
    function("logics_and",
             (SEXP (*)(const SEXP &lhs, const SEXP &rhs))&Logics::logics_and,
             "Builds the logical conjunction of two formulas");
    function("logics_or", (SEXP (*)(const SEXP &, const SEXP &))&Logics::logics_or,
             "Builds the logical disjunction of two formulas");
}
