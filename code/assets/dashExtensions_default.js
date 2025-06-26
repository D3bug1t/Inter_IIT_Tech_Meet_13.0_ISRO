window.dashExtensions = Object.assign({}, window.dashExtensions, {
    default: {
        function0: function(feature, context) {
            const {
                classes,
                colorscale,
                style,
                colorProp
            } = context.hideout
            style.fillColor = colorscale[feature.properties[colorProp]]; // set the fill color according to the class    
            return style;
        }
    }
});