(defparameter amt-of-drug '(1 2 3 4 5))
(defparameter reaction-time '(1 1 2 2 4))

(defparameter dist-from-fire-station
  '(3.4 1.8 4.6 2.3 3.1 5.5 .7 3.0 2.6 4.3 2.1 1.1 6.1 4.8 3.8))
(defparameter fire-damage
  '(26.2 17.8 31.3 23.1 27.5 36.0 14.1 22.3 19.6 31.3 24.0 17.3 43.2 36.4 26.1))

(defparameter maryj-opinions '(39 99 336 26))
(defparameter maryj-proportions '(.07 .18 .65 .1))

(defparameter lifestyles-of-marriage-p13.55
  (make-array '(5 5)
	      :initial-contents '((9 6 6 3 9)
				  (7 11 12 11 5)
				  (8 8 6 11 5)
				  (8 10 7 15 11)
				  (4 6 7 10 7))))

;;the rank is the position in sorted list
;;except for duplicates, each duplicate
;;adds .5 to the rank 
(defun wilcoxon-rank (l)
  (sort (copy-list l) #'<))

(defun large-sample-wilcoxon-rank-sum-test (t1 n1 n2)
  (/ (- t1 (/ (* n1 (+ n1 n2 1))
	      2))
     (sqrt (/ (* n1 n2 (+ n1 n2 1))
	      12))))

(defun lte-binomial-probability (n k p)
  (sum (mapcar (lambda (x) 
		 (exact-binomial-probability n x p))
	       (range k))))

(defun exact-binomial-probability (n k p)
  (* (expt p k)
     (expt (- 1 p) (- n k))
     (/ (fact n)
	(* (fact k) (fact (- n k))))))

(defun large-sample-sign-test-ts (s n)
  (/ (- (- s .5) (* .5 n))
     (* .5 (sqrt n))))

(defun chi^2-of-mat (mat)
  (loop
     for i from 0 to (1- (car (array-dimensions mat)))
     sum (chi^2-of-row mat i) into total
     finally (return total)))

(defun chi^2-of-row (mat row)
  (loop
     for i from 0 to (1- (cadr (array-dimensions mat)))
     sum (chi^2-of-cell mat row i) into row-total
     finally (return row-total)))

(defun chi^2-of-cell (mat row col)
  (/ (^2 (- (aref mat row col) (expected-cell mat row col)))
     (expected-cell mat row col)))

(defun expected-cell (mat row col)
  (/ (* (sum-row mat row) (sum-col mat col))
     (sum-mat mat)))
	
(defun sum-mat (mat)
  (loop
     for i from 0 to (1- (car (array-dimensions mat)))
     sum (sum-row mat i) into total
     finally (return total)))

(defun sum-row (mat row)
  (loop
     for i from 0 to (1- (cadr (array-dimensions mat)))
     sum (aref mat row i) into total
     finally (return total)))

(defun sum-col (mat col)
  (loop
     for i from 0 to (1- (car (array-dimensions mat)))
     sum (aref mat i col) into total
     finally (return total)))

(defun values-from-percentages (n percents)
  (mapcar #'* percents (make-list 5 :initial-element n)))

(defun chi^2-one-way (values expected-values)
  (sum (mapcar (lambda (n e) (/ (^2 (- n e))
				e))
	       values
	       expected-values)))

(defun simple-linear-regression (x y)
  (format t "~&y = Beta0 + Beta1(x)")
  (format t "~&Beta1 = ~v$" 3 (slr-beta1 x y))
  (format t "~&Beta0 = ~v$" 3 (slr-beta0 x y))
  (format t "~&y = ~V$ + ~V$(x)" 
	  3 (slr-beta0 x y)
	  3 (slr-beta1 x y))
  (format t "~&SSE = ~V$" 3 (slr-sse x y))
  (format t "~&MSE = ~V$" 3 (slr-mse x y)))

(defun slr-std-err-y-hat (x y specific-x)
  (sqrt (* (slr-mse x y)
	   (+ (/1 (length x))
	      (/ (^2 (- specific-x (mean x)))
		 (SSxx x))))))

(defun slr-std-err-prediction (x y specific-x)
  (sqrt (* (slr-mse x y)
	   (+ 1
	      (/1 (length x))
	      (/ (^2 (- specific-x (mean x)))
		 (SSxx x))))))

(defun slr-y-hat (x y specific-x)
  (+ (slr-beta0 x y)
     (* (slr-beta1 x y) specific-x)))

;;not done like the one below
(defun slr-prediction-interval (x y t-val specific-x)
  (- (slr-y-hat x y specific-x)
     (* t-val (slr-std-err-prediction x y specific-x))))


;;not done like the one below either
(defun slr-confidence-interval-at-x (x y t-val specific-x)
  (- (slr-y-hat x y specific-x)
     (* t-val (slr-std-err-y-hat x y specific-x))))



(defun slr-r^2 (x y)
  (- 1 (/ (slr-sse x y)
	  (ssyy y))))

(defun slr-r (x y)
  (/ (SSxy x y)
     (sqrt (* (SSxx x) (SSyy y)))))

(defun slr-beta1-confidence-interval (x y t-val)
  (format t "~&~V$ +- ~$ * ~V$"
	  3 (slr-beta1 x y)
	  t-val
	  3 (slr-std-err-beta1 x y))
  (format t "~&(~V$, ~V$)" 
	  3 (- (slr-beta1 x y) 
	     (* t-val (slr-std-err-beta1 x y)))
	  3 (+ (slr-beta1 x y) 
	     (* t-val (slr-std-err-beta1 x y)))))

;;estimated standard error of beta1-hat
(defun slr-std-err-beta1 (x y)
  (/ (sqrt (slr-mse x y))
     (sqrt (SSxx x))))

(defun slr-mse (x y)
  (/ (slr-sse x y) (- (length x) 2)))

(defun slr-sse (x y)
  (sum (mapcar (lambda (xi yi)
		 (^2 (- yi (+ (slr-beta0 x y)
			      (* xi (slr-beta1 x y)))))) 
	       x y)))

(defun slr-beta1 (x y)
  (/ (SSxy x y) (SSxx x)))

(defun slr-beta0 (x y)
  (- (mean y) (* (mean x) (slr-beta1 x y))))

(defun SSyy (y)
  (SSxy y y))

(defun SSxx (x)
  (SSxy x x))

(defun SSxy (x y)
  (sum (mapcar (lambda (xi yi) 
		 (* (- xi (mean x))
		    (- yi (mean y))))
	       x y)))

(defun mean (x)
  (/ (reduce #'+ x) (length x)))

(defun ^2 (x) (expt x 2))
(defun ^3 (x) (expt x 3))

(defun /1 (x) (/ 1 x))

(defun fact (x)
  (cond ((< x 1)
	 1)
	(t
	 (* x (fact (1- x))))))

(defun range (x)
  (loop for i upto x collect i))

(defun sum (x) (reduce #'+ x))

(defun interval-printer (left right)
  (format t "~&~V$ +- ~V$"
	  3 left
	  3 right)
  (format t "~&(~V$, ~V$)" 
	  3 (- left
	       right)
	  3 (+ left
	       right)))