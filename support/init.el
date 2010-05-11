;; This adds additional extensions which indicate files normally
;; handled by cc-mode.
(setq auto-mode-alist
      (append '(("\\.C$"  . c++-mode)
		("\\.cc$" . c++-mode)
		("\\.hh$" . c++-mode)
		("\\.h$"  . c++-mode))
	      auto-mode-alist))

;; Set stroustrup as the default style for C/C++ code
(setq c-default-style "stroustrup")

;; Set up C++ mode hook
(defun my-c++-mode-hook ()
  ;; Tell cc-mode not to check for old-style (K&R) function declarations.
  ;; This speeds up indenting a lot.
  (setq c-recognize-knr-p nil)

  ;; switch/case:  make each case line indent from switch
  (c-set-offset 'case-label '+)

  ;; Automatically indent after a newline (like vi)
  (local-set-key '[(return)] 'newline-and-indent)

  ;; Tab sanity: match the C indent, but don't insert new tabs (use spaces)
  (setq tab-width 4)
  (setq indent-tabs-mode nil))

(add-hook 'c++-mode-hook 'my-c++-mode-hook)

