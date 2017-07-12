var path = require('path');
var webpack = require('webpack');

module.exports = {
  entry: {
    app: './src/main.js',
    vendor: ['d3']
  },
  plugins: [
    new webpack.optimize.CommonsChunkPlugin('vendor', 'dst/vendor.bundle.js'),
    new webpack.optimize.UglifyJsPlugin({
      compress: { warnings: false }
    }),
    new webpack.NoErrorsPlugin(),
  ],
  output: {
    path: __dirname,
    filename: 'dst/bundle.js'
  },
  module: {
    preLoaders: [
      {
        test: /\.js$/,
        loader: "eslint-loader",
        exclude: /node_modules/
      }
    ],
    loaders: [
      {
        loader: 'babel-loader',
        exclude: /node_modules/,
        query: {
          presets: ['es2015']
        }
      }
    ]
  },
};
