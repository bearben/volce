(set-logic QF_LIA)
(set-info :smt-lib-version 2.0)
(set-info :status sat)
(declare-fun c0 () Int)
(declare-fun c1 () Int)
(declare-fun c2 () Int)
(assert (or (= c0 32) (= c0 9) (= c0 10)))
(assert (not (or (= c1 32) (= c1 9) (= c1 10))))
(assert (not (and (distinct c1 46) (or (< c1 48) (> c1 57)))))
(assert (not (and (>= c2 48) (<= c2 57))))
(assert (not (= c2 46)))
(check-sat)
(exit)
