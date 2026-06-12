====================================================================================
Observations and lessons learned so far on scaling and performance in PyPSA
====================================================================================

1. Solver-side scaling and model-side scaling
=============================================

1.1 Solver-side scaling
-------------------

Implemented internally by the solver. Its purpose is to improve the numerical conditioning of
the optimization problem before or during the solution process.

Matrix scaling techniques
^^^^^^^^^^^^^^^^^^^^^^^^^

Several matrix scaling techniques are used in optimization solvers to improve the numerical conditioning of the constraint matrix while preserving the original optimization problem.

Given the linear system:

.. math:: Ax = b,

where

* :math:`A \in \mathbb{R}^{m \times n}` is the **constraint matrix** containing the coefficients of the optimization problem,
* :math:`x \in \mathbb{R}^{n}` is the vector of **decision variables**,
* :math:`b \in \mathbb{R}^{m}` is the **right-hand-side (RHS) vector**,

diagonal row and column scaling matrices can be introduced:

.. math:: R =\mathrm{diag}(r_1, r_2, \ldots, r_m),

.. math:: C = \mathrm{diag}(c_1, c_2, \ldots, c_n),

where

* :math:`R \in \mathbb{R}^{m \times m}` is the **row scaling matrix**, whose diagonal entries rescale the constraints (rows) of the optimization problem,
* :math:`C \in \mathbb{R}^{n \times n}` is the **column scaling matrix**, whose diagonal entries rescale the decision variables (columns) of the optimization problem.

The scaled decision variables :math:`y` are defined through the transformation

.. math:: x = Cy,

where :math:`y \in \mathbb{R}^{n}` denotes the vector of **scaled decision variables**.

Substituting this transformation into the original system and applying row scaling yields the scaled problem:

.. math:: RACy = Rb.

Equilibration scaling
^^^^^^^^^^^^^^^^^^^^^

Equilibration scaling aims to balance the magnitudes of coefficients across rows and columns of the constraint matrix. The scaling factors are typically computed iteratively such that rows and columns exhibit similar norms. A common approach relies on infinity norms:

.. math:: # R_{ii}\frac{1}{\max_j |a_{ij}|},

.. math:: # C_{jj}\frac{1}{\max_i |(RA)_{ij}|}.

The scaling process may be repeated until the matrix becomes sufficiently balanced.

Max value scaling
^^^^^^^^^^^^^^^^^

Max value scaling normalizes each row and/or column using its largest absolute coefficient. For a row/column :math:`i`, the scaling factors are defined as:

.. math:: # R_{ii}\frac{1}{\max_j |a_{ij}|},
.. math:: # C_{jj}\frac{1}{\max_i |a_{ij}|}.

As a result, the largest coefficient in each scaled row or column becomes approximately equal to one.

Geometric scaling
^^^^^^^^^^^^^^^^^

Geometric scaling determines scaling factors using both the smallest and largest non-zero coefficients. For a row :math:`i`, the scaling factor can be expressed as

.. math:: # R_{ii}\frac{1}{\sqrt{\min_{j:a_{ij}\neq0}|a_{ij}|\max_{j:a_{ij}\neq0}|a_{ij}|}}.

An analogous expression can be defined for columns.

By relying on the geometric mean of the extreme coefficient magnitudes, this approach seeks to reduce the spread of coefficients while avoiding an excessive focus on the largest entries.

Objective scaling
^^^^^^^^^^^^^^^^^

Objective scaling modifies the objective coefficients :math:`\min c^\top x` by replacing :math:`c` with :math:`\tilde{c} = \alpha c`, where :math:`\alpha > 0` is a scaling factor.

Since the objective function is multiplied by a positive constant, the optimal solution remains unchanged. Some solvers, such as HiGHS, implement this using powers of two, :math:`\alpha = 2^k` to avoid introducing additional floating-point round-off errors.

Bound scaling
^^^^^^^^^^^^^

Bound scaling modifies variable bounds, :math:`l \leq x \leq u`, according to:

.. math:: \tilde{l} = \beta l,
.. math:: \tilde{u} = \beta u,

where :math:`\beta > 0` is a scaling factor.

This operation preserves the feasible region after the corresponding variable transformation and may improve numerical behaviour when variable bounds differ substantially in magnitude.

Remarks
^^^^^^^^

These techniques improve the numerical conditioning of the optimization problem but do not alter its underlying mathematical structure. In particular, they cannot compensate for poor model formulation, such as excessively large big-:math:`M` constants, unrealistic variable bounds, or inconsistent modelling units.

The exact implementation of these scaling techniques varies across solvers, and the equations above represent their general principles rather than solver-specific algorithms.

Solvers such as **Gurobi** automatically select an internal scaling strategy to improve the numerical conditioning of the optimization problem, typically relying on approaches related to equilibration or geometric scaling.

**HiGHS** also performs internal matrix scaling by default. Its default approach is equilibration scaling, which attempts to balance the magnitude of coefficients across the constraint matrix. Alternative scaling strategies can be selected by the user, including disabling scaling entirely or adopting max value scaling.

1.2 Model-side scaling
------------------

This happens during model formulation.

The goal is to **construct an optimization problem that is numerically
well-behaved from the outset**.

Examples include:

-  using consistent units,
-  avoiding unnecessarily large or small coefficients,
-  avoiding excessively loose bounds,
-  removing unnecessary constants,
-  rescaling constraints when appropriate.

This is the area where **PyPSA could potentially provide native
support**.

During the OET webinar, both HiGHS and Gurobi developers emphasized that
model scaling is extremely important.

The recommended target is generally to have **coefficients around O(1)** rather than having coefficients spanning many orders of magnitude.

On their hand, **solver developers are already working extensively on internal scaling strategies**. Therefore, improvements within PyPSA should probably focus on helping users formulate better-scaled models, rather than attempting to outperform solver-internal scaling mechanisms.

Using a limited set of PyPSA-Eur LP power-only models (``pypsa-de-elec`` in different spatial resolutions and 1h temporal resolution), all the solver-side scaling strategies implemented in Gurobi fail in solving the following warning:

::

   Warning: Model contains large bounds

while the reported coefficient statistics were, e.g.:

::

   Matrix range     [1e-02, 8e+01]
   Objective range  [9e-03, 3e+05]
   Bounds range     [1e+10, 1e+10]
   RHS range        [1e+00, 3e+05]

The electricity-only examples lead to the important conclusion that **the constraint matrix itself was already reasonably well scaled**. Therefore, solver-side scaling is likely already doing its job adequately.

For instance, increasing the model size from 10 to 50 nodes while maintaining hourly resolution yielded matrix coefficient ranges approximately between :math:`10^{-2}` and :math:`10^{2}`, despite a substantial increase in the number of variables and constraints.

However, the persistent warning concerned large variable bounds:

::

   Bounds range [1e+10, 1e+10]

Additional sector-coupled examples (``pypsa-de-sec``) confirmed that matrix coefficient ranges remained relatively well behaved, approximately between :math:`10^{-3}` and :math:`10^{2}`, while introducing substantially larger RHS values, reaching up to :math:`10^{9}`.

These observations suggest that the main numerical concerns are unlikely to originate from poor constraint matrix scaling. Instead, they appear to be related to model formulation, particularly extremely loose variable bounds and, in sector-coupled models, very large RHS values.

Potential sources include:

-  variables with effectively infinite bounds;
-  large default upper bounds replacing infinity;
-  big-M (large coefficients intentionally introduced in the formulation to model conditional relationships) formulations;
-  unconstrained investment capacities.

In conclusion, only the model formulation has access to the physical meaning needed to determine whether a bound of ``1e10`` or an RHS value of :math:`10^9` is reasonable.

2. Scaling diagnostics within linopy
====================================

Optimization solvers already report aggregate statistics describing the numerical properties of the optimization problem. Nevertheless, these diagnostics are often difficult to interpret from a modelling perspective.

PyPSA-native diagnostics could complement solver outputs by linking numerical issues back to the underlying model formulation. Rather than simply reporting that large coefficient ranges or variable bounds exist, such diagnostics could identify the specific variables, constraints, and components responsible for these characteristics.

Potential diagnostics could include:

- identification of variable groups associated with large bounds,
- identification of constraints with large internal coefficient ratios,
- summaries of objective contributions spanning multiple orders of magnitude,
- warnings highlighting modelling patterns commonly associated with numerical difficulties.

This additional layer of interpretation could help users diagnose and address scaling issues directly at the modelling level.

3. AnyMod approach
==================

The AnyMod implementation already performs scaling directly within the model.

The procedure essentially:

#. computes coefficient ranges for each constraint;
#. identifies constraints outside a desired range;
#. rescales the entire constraint expression by a scalar factor.

This corresponds primarily to **constraint row scaling**.

Importantly:

-  it preserves feasibility,
-  it preserves optimal solutions,
-  it changes only the numerical representation.

However, it does **not** address:

-  poor variable bounds,
-  column scaling,
-  problematic objective scaling.

4. Proposed pathway
===================

A reasonable development pathway for PyPSA could be:

Step 1: Diagnostics
-------------------

Via Linopy, provide users with diagnostics that link solver-reported numerical characteristics back to the underlying PyPSA formulation.

Potential outputs could include:

- variables contributing to excessively large bounds,
- constraints responsible for large RHS values,
- identification of constraints with extreme coefficient ratios,
- summaries of objective terms spanning multiple orders of magnitude,
- identification of the model components responsible for these characteristics.

Step 2: Safe formulation improvements
-------------------------------------

In PyPSA:

- avoid artificial infinities where possible,
- introduce tighter, physically meaningful bounds where appropriate,
- investigate whether additional formulation improvements preserving the mathematical problem can be identified.

PyPSA already includes (PR [#1509](https://github.com/PyPSA/PyPSA/pull/1509)) as an example of model-side numerical improvement. Before that, constant objective contributions associated with pre-existing capacities unnecessarily increased the magnitude of the objective function without affecting the optimization problem itself.

In particular:
#. it introduced the ``include_objective_constant`` parameter to both ``n.optimize()`` and ``n.optimize.create_model()``. This parameter controls whether the objective constant, corresponding to the capital costs associated with existing infrastructure on extendable assets, is explicitly included in the optimization objective.
#. When ``include_objective_constant=None`` (the current default), the objective constant is included in the optimization model. Since this contribution is independent of the decision variables, it does not affect the optimal solution but increases the magnitude of the objective function. Setting ``include_objective_constant=False`` excludes this constant contribution from the optimization objective, thereby improving the numerical representation of the LP without altering the feasible set or optimal solution.
#. The objective constant is still computed and stored on the network object (``n.objective_constant``) regardless of the value of ``include_objective_constant``. Therefore, this modification affects only the formulation of the optimization problem and not the reported total system costs.

All in all, the implemented solution separated these constant terms from the optimization problem by treating them as an objective constant rather than including them directly in the optimization objective.

This modification preserved the feasible set and optimal solution while reducing unnecessary numerical burden associated with large objective values.

Step 3: Optional row scaling
----------------------------

Potentially inspired by AnyMod, but clearly documented and disabled by default.

Given the evidence collected so far, however, optional row scaling should be regarded as an exploratory improvement rather than a primary target, since the tested PyPSA models already exhibit relatively well-scaled constraint matrices.

5. On Benders decomposition
===========================

Göke et al. (*Stabilized Benders decomposition for energy planning under climate uncertainty*)[https://www.sciencedirect.com/science/article/pii/S0377221724000353?via%3Dihub] propose the application of Benders decomposition to large-scale energy system planning problems under uncertainty.

Benders decomposition addresses computational complexity by decomposing a large optimization problem into multiple smaller problems that are solved iteratively. Typically, a **master problem** determines a subset of decision variables (e.g. investment decisions), while one or more **subproblems** evaluate the feasibility and optimality of these decisions under different operational conditions or scenarios. Information from the subproblems is then communicated back to the master problem through additional constraints, commonly referred to as *Benders cuts*.

Unlike scaling techniques, **Benders decomposition does not aim to improve the numerical conditioning of the optimization problem**. Instead, it **seeks to improve computational tractability by exploiting the underlying structure of the problem and reducing the effective size of the optimization problems solved at each iteration**.

The stabilized variant proposed by Göke et al. further improves the convergence behaviour of the decomposition algorithm by reducing oscillations in the master problem between successive iterations. In this context, *stabilization* refers to improved algorithmic convergence rather than enhanced numerical conditioning.

Consequently, Benders decomposition should be regarded as a separate class of performance-improvement techniques. Although decomposition approaches may indirectly alleviate some numerical difficulties by operating on smaller subproblems, their primary objective is not numerical robustness but computational efficiency.

In the context of PyPSA, decomposition approaches would therefore represent an alternative strategy for improving solver performance, complementary to both solver-side scaling techniques (already applied consistently if using Gurobi or comparable commercial solvers) and model-side numerical improvements.

5. Conclusion
=============

So far, solver-side scaling is already handled extensively by solver developers.

Consequently, the most valuable contribution PyPSA could make is likely in the area of model-side support, including:

-  numerical diagnostics,
-  identification of problematic bounds,
-  formulation guidance,
-  optional safe rescaling procedures.

In particular, the examples explored indicate that large bounds may currently represent a more significant source of numerical difficulty than the scaling of the constraint matrix itself.

6. Open questions
=================

- Investigate the origin of large default bounds (e.g. :math:`10^{10}`) and assess whether tighter, physically meaningful bounds can be introduced.
- Investigate the origin of large RHS values observed in sector-coupled models.
- Explore PyPSA-native diagnostics that attribute solver-reported numerical issues to specific model components and formulations.
- Evaluate whether additional model-side numerical improvements that preserve the mathematical formulation can be identified.
- Investigate whether optional constraint row scaling, inspired by the approach adopted in AnyMod, provides measurable benefits for PyPSA models.
- Keep decomposition and problem-size reduction approaches (e.g. Benders decomposition or Issue #1589) conceptually separate from scaling-related techniques.

