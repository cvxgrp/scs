.. _relaxation:

Relaxation
==========

The (unscaled, see :ref:`scaling`) SCS update equations are:

.. math::
  \begin{align}
  \tilde u^{k+1} &= (I + \mathcal{Q})^{-1} w^k \\
  u^{k+1} &= (I + N_{\mathcal{C}_+})^{-1} (2 \tilde u^{k+1} - w^k) \\
  w^{k+1} &= w^k + \alpha (u^{k+1} - \tilde u^{k+1}) \\
  \end{align}

where :math:`\alpha \in (0,2)` is the *relaxation* parameter. Vanilla
Douglas-Rachford corresponds to setting :math:`\alpha = 1`. If :math:`\alpha <
1` it is referred to as under-relaxation, if :math:`\alpha > 1` it is
over-relaxation.  Typically values of :math:`\alpha \approx 1.5` work well.  It
is controlled by the :code:`alpha` :ref:`setting <settings>`.  If :math:`\alpha
= 2` the method reduces to Peaceman-Rachford splitting, which is not guaranteed
to converge in general (but will if, for example, the operators in the problem
are both maximal strongly monotone).

Thankfully, there is no interaction between :math:`\alpha` and the scaling
described in :ref:`scaling`, and they can be combined immediately.
