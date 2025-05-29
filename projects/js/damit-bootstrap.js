"use strict";

/** Enable popovers everywhere. */
$(function() {
	// If delay is set as attribute, don't override it.
	$('[data-toggle="popover"][data-delay]').popover({ container: "body" });
	$('[data-toggle="popover"]').popover({ delay: 300, container: "body" });
});

/** Close popover when clicking outside. (https://jsfiddle.net/mattdlockyer/C5GBU/2/) */
$("body").on("click", function(e) {
	$('[data-toggle="popover"]').each(function() {
		//the 'is' for buttons that trigger popups
		//the 'has' for icons within a button that triggers a popup
		if (
			!$(this).is(e.target) &&
			$(this).has(e.target).length === 0 &&
			$(".popover").has(e.target).length === 0
		) {
			$(this).popover("hide");
		}
	});
});

/** Put full text into title when .text-truncate is used and element is overflowing. */
$(function() {
	$(".text-truncate").each(function() {
		let elem = $(this);
		if (
			elem.prop("scrollHeight") > elem.height() ||
			elem.prop("scrollWidth") > elem.width()
		) {
			elem.attr("title", elem.text());
		}
	});
});
