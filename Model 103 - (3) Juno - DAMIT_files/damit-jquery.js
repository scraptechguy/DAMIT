"use strict";

/** Center an element over another element or parent. */
jQuery.fn.extend({
	centerOver: function (other) {
		other = (typeof other !== "undefined") ? other : $(this).parent();
		let otherTop = $(other).position().top;
		let otherLeft = $(other).position().left;
		let otherHeight = $(other).height();
		let otherWidth = $(other).width();
		let top = otherTop + (otherHeight - $(this).height())/2;
		let left = otherLeft + (otherWidth - $(this).width())/2;
		this.css("position", "absolute");
		this.css("top", top + "px");
		this.css("left", left + "px");
		return this;
	}
});
